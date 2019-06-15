#define TMB_LIB_INIT R_init_spatq
#include <TMB.hpp>
#include "../inst/include/spatq.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_VECTOR(catch_obs);
  DATA_MATRIX(X_n);
  DATA_MATRIX(X_w);
  DATA_STRUCT(spde, spde_t);
  DATA_SPARSE_MATRIX(A);

  PARAMETER_VECTOR(beta_n);
  PARAMETER_VECTOR(beta_w);
  PARAMETER_VECTOR(log_kappa);
  PARAMETER_VECTOR(log_tau);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(spat_n);
  PARAMETER_VECTOR(spat_w);

  // Get number of observations for later
  int N_obs = catch_obs.size();

  // Initialize joint negative log-likelihood accumulator
  vector<Type> jnll(3);
  jnll.setZero();

  vector<Type> kappa = exp(log_kappa);
  vector<Type> tau = exp(log_tau);
  Type sigma = exp(log_sigma);

  // Project spatial effects from mesh nodes to observation locations
  vector<Type> omega_n(N_obs);
  vector<Type> omega_w(N_obs);
  omega_n = (A * spat_n) / tau(0);
  omega_w = (A * spat_w) / tau(1);

  // Get density of spatial random effects, repeated for each process
  SparseMatrix<Type> Q_n = Q_spde(spde, kappa(0));
  SparseMatrix<Type> Q_w = Q_spde(spde, kappa(1));
  jnll(1) += GMRF(Q_n)(spat_n);
  jnll(2) += GMRF(Q_w)(spat_w);

  // Get n and w for each observation
  vector<Type> n(N_obs);
  vector<Type> w(N_obs);
  n = X_n * beta_n + omega_n;
  w = X_n * beta_n + omega_n;

  // Convert to p and r
  vector<Type> p(N_obs);
  vector<Type> r(N_obs);
  vector<Type> pr_tmp(2);
  for (int i = 0; i < N_obs; i++) {
    pr_tmp = poislink(n(i), r(i));
    p(i) = pr_tmp(0);
    r(i) = pr_tmp(1);
  }

  // Calculate likelihood
  for (int i = 0; i < N_obs; i++) {
    if (catch_obs(i) == 0) {
      jnll(0) -= log(Type(1.0) - p(i));
    } else {
      jnll(0) -= log(p(i)) +
        dlnorm(catch_obs(i), r(i) - sigma * sigma / Type(2.0), sigma, true);
    }
  }

  REPORT(jnll);

  return sum(jnll);
}

