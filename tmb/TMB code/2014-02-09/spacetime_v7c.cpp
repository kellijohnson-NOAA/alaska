// Space time 
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(n_data);         // Total number of observations
  DATA_VECTOR(Y);       	// Count data
  DATA_FACTOR(NAind);		// 1 = Y is NA, 0 = is not NA
  DATA_FACTOR(NegBin);      // Use negative binomial? 0:No; 1:Yes
  DATA_INTEGER(n_stations)	// Number of stations 
  DATA_FACTOR(meshidxloc);	// Pointers into random effects vector x
  DATA_INTEGER(n_years)          // Number of years  
  DATA_INTEGER(n_p)          	// number of columns in covariate matrix X
  DATA_MATRIX(X);		// Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  PARAMETER_VECTOR(alpha);   // Mean of Gompertz-drift field
  PARAMETER(phi);            // Offset of beginning from equilibrium
  PARAMETER(log_tau_E);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)
  PARAMETER_VECTOR(ln_VarInfl);   // Overdispersion parameters

  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in carrying capacity

  using namespace density;
  int i,j;
  Type g = 0;
  
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  using namespace density;
  g += SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);
  g += GMRF(Q)(Omega_input);

  // Likelihood contribution from observations
  vector<Type> eta(n_data); 
  vector<Type> nu(n_data);
  vector<Type> mean_abundance(n_years);
  matrix<Type> Dji(n_stations,n_years);
  matrix<Type> Epsilon(n_stations,n_years);
  vector<Type> Omega(n_stations);
  vector<Type> Equil(n_stations);

  eta = X*alpha.matrix();
  int ii = 0;
  for (int j=0;j<n_stations;j++){
    Omega(j) = Omega_input(j) / exp(log_tau_O);
    Equil(j) = eta(ii) + Omega(j) / (1-rho);
    ii++;
  }
  
  ii = 0;
  for (int i=0;i<n_years;i++)
  {
    //mean_ln_abundance(i) = 0;
    mean_abundance(i) = 0;
    for (int j=0;j<n_stations;j++)
    { 
      Epsilon(j,i) = Epsilon_input(j,i) / exp(log_tau_E);
      nu[ii] = phi * pow(rho, i);
      Dji(j,i) = nu[ii] + Epsilon_input(meshidxloc(j),i)/exp(log_tau_E) + ( eta(ii) + Omega_input(meshidxloc(j))/exp(log_tau_O) )/(1-rho);
      mean_abundance[i] = mean_abundance[i] + exp( Dji(j,i) );      
      if(!NAind(ii)){                
        Type mean_y = exp( Dji(j,i) );
        Type var_y = mean_y*(1.0+exp(ln_VarInfl(0)))+pow(mean_y,2.0)*exp(ln_VarInfl(1));
        if(!NegBin(0)){ g -= dpois( Y[ii], mean_y, 1 ); }
        if(NegBin(0)){ g -= dnbinom( Y[ii], mean_y, var_y, 1 ); }
      }
      ii++;      
    }
    mean_abundance[i] = mean_abundance[i] / n_stations;      
  }

  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  // Fields
  REPORT( Dji );
  REPORT( Epsilon );
  REPORT( Omega );
  REPORT( Equil );
  // Total abundance
  ADREPORT( log(mean_abundance) );
  ADREPORT(mean_abundance);
  
  return g;
}
