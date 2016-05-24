// Space time
#include <TMB.hpp>

/* Detecting NAs */
template <class Type>
bool isNA(Type x) {
  return R_IsNA(asDouble(x));
}

/* Parameter transform */
template <class Type>
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(c_i);       	   // Count data
  DATA_VECTOR(year);           // values of 1:24 for every c_i
  DATA_VECTOR(station_map);    // Points into RE vector for every pred
                               // should only be used with Dji
  DATA_VECTOR(station_unique); // Points to used portions of RE vector
                               // should only be used with "_input"
  DATA_INTEGER(n_x)            // Number of stations
  DATA_INTEGER(n_years)        // Number of years including init

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  PARAMETER_VECTOR(alpha);   // Mean of Gompertz-drift field
  PARAMETER(phi);            // Offset of beginning from equilibrium
  PARAMETER(log_tau_E);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(log_sigma);
  PARAMETER(rho);            // Autocorrelation (i.e. density dependence)

  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in carrying capacity

  using namespace density;
  int i, ii, j, station_use, year_use;
  vector<Type> g(3);

  Type kappa2 = exp(2.0 * log_kappa);
  Type kappa4 = kappa2 * kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp(log_kappa);
  Type SigmaE = 1 / sqrt(4 * pi * exp(2 * log_tau_E) * exp(2 * log_kappa));
  Type SigmaO = 1 / sqrt(4 * pi * exp(2 * log_tau_O) * exp(2 * log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  using namespace density;
  // Compute the random effects matrix for process error
  // initial year is included in the [0] index
  g(0) = 0;
  for (int i = 0; i < n_years; i++){
    g(0) += GMRF(Q)(Epsilon_input.col(i));
  }
  // Compute the random effects vector for the spatial field
  g(1) = GMRF(Q)(Omega_input);

  // Likelihood contribution from observations
  vector<Type> mean_abundance(n_years);
  matrix<Type> Dji(n_x, n_years);
  // May want to remove lnB_stations as it is just the transpose
  // of Dji
  matrix<Type> lnB_stations(n_years, n_x);
  matrix<Type> Epsilon(n_x, n_years);
  vector<Type> Omega(n_x);
  vector<Type> Equil(n_x);

  // calculate theoretical equilibrium for each
  // node that is occupied
  ii = 0;
  for (int j = 0; j < n_x; j++){
    int station_use = CppAD::Integer(station_unique(j));
    Omega(j) = Omega_input(station_use) / exp(log_tau_O);
    Equil(j) = alpha(0) + Omega(j) / (1 - rho);
  }

  // Initial year
  for (int j = 0; j < n_x; j++){
    int station_use = CppAD::Integer(station_unique(j));
    Epsilon(j, 0) = Epsilon_input(station_use, 0) / exp(log_tau_E);
    Dji(j, 0) = phi + (alpha(0) + Omega_input(station_use) / exp(log_tau_O)) / (1 - rho) + Epsilon_input(station_use, 0) / exp(log_tau_E);
  }
  // Recursive formula for subsequent years
  for (int i = 1; i < n_years; i++){
    for (int j = 0; j < n_x; j++){
      int station_use = CppAD::Integer(station_unique(j));
      Epsilon(j, i) = Epsilon_input(station_use, i) / exp(log_tau_E);
      Dji(j, i) = alpha(0) + Omega_input(station_use) / exp(log_tau_O) + Dji(j, i - 1) * rho + Epsilon_input(station_use, i) / exp(log_tau_E);
    }
  }

  // calculate the mean abundance over all occupied stations
  // in each given year, including predicted init year
  // mean abundance is on the same scale as the data
  g(2) = 0;
  for (int i = 0; i < n_years; i++){
    mean_abundance(i) = 0;
    for (int j = 0; j < n_x; j++){
      lnB_stations(i, j) = Dji(j, i);
      mean_abundance[i] = mean_abundance[i] + exp(Dji(j, i));
    }
    mean_abundance[i] = mean_abundance[i] / n_x;
  }

  // Probability of data
  // year ranges from 1:24 because year 0 has no observed data
  for(int ii = 0; ii < c_i.size(); ii++){
    int station_use = CppAD::Integer(station_map(ii));
    int year_use = CppAD::Integer(year(ii));
    Type mean_y = Dji(station_use, year_use);
    Type sigma_y = exp(log_sigma);
    g(2) -= log(1 / (c_i[ii] * sigma_y * sqrt(2.0*M_PI)) * exp(-square(log(c_i[ii]) - mean_y) / (2 * square(sigma_y))));
    //g(2) += 0.5*(log(2.0*M_PI*sigma_y) + square(log(c_i[ii]) - mean_y)/sigma_y);
    //g(2) -= dnorm(c_i[ii], mean_y, sigma_y, 1);
  }

  // Spatial field summaries
  REPORT(Range);
  REPORT(SigmaE);
  REPORT(SigmaO);
  // Fields
  REPORT(Epsilon);
  REPORT(Omega);
  REPORT(Equil);
  REPORT(g);
  REPORT(log_sigma);
  // Total abundance
  ADREPORT(log(mean_abundance));
  ADREPORT(mean_abundance);
  ADREPORT(Dji);

  Type g_sum = sum(g);
  return g_sum;
}
