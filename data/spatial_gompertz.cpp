// Space time
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type log_mean, Type log_sd, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), log_mean, exp(log_sd), false) / x;
  if(give_log==true) Return = dnorm( log(x), log_mean, exp(log_sd), true) - log(x);
  return Return;
}

template<class Type>
Type d_poisson_lognormal(Type x, Type log_mean, Type log_sd, Type log_clustersize, int give_log=false){
  Type Return;
  Type log_notencounterprob = -1 * exp(log_mean) / exp(log_clustersize);
  Type encounterprob = 1 - exp( log_notencounterprob );
  if( x==0 ){
    Return = log_notencounterprob;
  }else{
    Return = log(encounterprob) + dlognorm( x, log_mean-log(encounterprob), log_sd, true );
  }
  if( give_log==true){ return Return; }else{ return exp(Return); }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_IVECTOR( Options_vec );
  // Slot 0:  Observation model (0=Poisson;  1=Poisson-lognormal)

  // Indices
  DATA_INTEGER( n_i );         // Total number of observations
  DATA_INTEGER( n_x );         // Number of vertices in SPDE mesh
  DATA_INTEGER( n_t );         // Number of years
  DATA_INTEGER( n_p );         // Number of columns in covariate matrix X

  // Data
  DATA_IVECTOR( x_s );	      // Association of each station with a given vertex in SPDE mesh
  DATA_VECTOR( c_i );       	// Count data
  DATA_IVECTOR( s_i );        // Station for each sample 0:(# of used unique stations - 1)
  DATA_IVECTOR( t_i );        // Time for each sample
  DATA_MATRIX( X_xp );		    // Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(alpha);   // Mean of Gompertz-drift field
  PARAMETER(phi);            // Offset of beginning from equilibrium
  PARAMETER(log_tau_E);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);            // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(theta_z); // Parameters governing measurement error

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in carrying capacity

  // objective function -- joint negative log-likelihood
  using namespace density;
  Type jnll = 0;
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Spatial parameters
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Objects for derived values
  vector<Type> eta_x(n_x);
  array<Type> log_Dpred_xt(n_x, n_t); //todo: determine what this should be for
  vector<Type> Omega_x(n_x);
  vector<Type> Equil_x(n_x);
  matrix<Type> Epsilon_xt(n_x, n_t);

  // Probability of Gaussian-Markov random fields (GMRFs)
  jnll_comp(0) += GMRF(Q)(Omega_input);
  jnll_comp(1) = SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);

  // Transform GMRFs
  eta_x = X_xp * alpha.matrix();
  for(int x=0; x<n_x; x++){
    Omega_x(x) = Omega_input(x) / exp(log_tau_O);
    Equil_x(x) = ( eta_x(x) + Omega_x(x) ) / (1-rho);
    for( int t=0; t<n_t; t++){
      Epsilon_xt(x,t) = Epsilon_input(x,t) / exp(log_tau_E);
    }
  }

  // Likelihood contribution from observations
  vector<Type> log_chat_i(n_i);
  vector<Type> jnll_i(n_i);
  jnll_i.setZero();
  vector<Type> log_notencounterprob_i(n_i);
  vector<Type> encounterprob_i(n_i);
  vector<Type> tmp_i(n_i);
  for (int i=0; i<n_i; i++){
    // t_i(i) is actually the timestep - 1 b/c indexing starts at zero
    // rho^0 == 1, and thus the population starts at phi
    log_chat_i(i) = phi*pow(rho,t_i(i)) + Epsilon_xt(x_s(s_i(i)),t_i(i)) + (eta_x(x_s(s_i(i))) + Omega_x(x_s(s_i(i))) ) / (1-rho);
    if( !isNA(c_i(i)) ){
      if(Options_vec(0)==0) jnll_i(i) -= dpois( c_i(i), exp(log_chat_i(i)), true );
      if(Options_vec(0)==1) jnll_i(i) -= d_poisson_lognormal( c_i(i), log_chat_i(i), theta_z(0), theta_z(1), true );
    }
  }
  jnll_comp(2) = jnll_i.sum();
  jnll = jnll_comp.sum();

  // Diagnostics
  REPORT( jnll_comp );
  REPORT( jnll );
  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( rho );
  ADREPORT( Range );
  ADREPORT( SigmaE );
  ADREPORT( SigmaO );
  // Fields
  REPORT( Epsilon_xt );
  REPORT( Omega_x );
  REPORT( Equil_x );
  // Diagnostics
  REPORT( log_chat_i );
  REPORT( jnll_i );
  REPORT( theta_z );
  REPORT(x_s);
  REPORT(c_i);
  REPORT(s_i);
  REPORT(t_i);
  REPORT(alpha);
  REPORT(phi);
  REPORT(log_tau_E);
  REPORT(log_tau_O);
  REPORT(log_kappa);
  REPORT(eta_x);
  REPORT(log_Dpred_xt);
  REPORT(log_notencounterprob_i);
  REPORT(encounterprob_i);

  return jnll;
}
