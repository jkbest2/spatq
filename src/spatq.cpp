#define TMB_LIB_INIT R_init_spatq
#include <TMB.hpp>
#include "../inst/include/spatq.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_VECTOR(catch_obs);         // observed catches
  DATA_MATRIX(X_n);               // fixed effect design matrices
  DATA_MATRIX(X_w);
  DATA_STRUCT(spde, spde_t);      // FEM matrices for spatial effects
  DATA_SPARSE_MATRIX(A);          // Projection matrix from mesh to obs

  PARAMETER_VECTOR(beta_n);       // mean log group density
  PARAMETER_VECTOR(beta_w);       // mean log weight per group
  PARAMETER_VECTOR(log_kappa);    // kappa parameter for each spatial field
  PARAMETER_VECTOR(log_tau);      // tau parameter for each spatial field
  PARAMETER(log_sigma);           // log catch variation
  PARAMETER_VECTOR(spat_n);       // spatial random effects on log group density
  PARAMETER_VECTOR(spat_w);       // spatial on log group weight

  // Get number of observations for later
  int N_obs = catch_obs.size();

  // Initialize joint negative log-likelihood accumulator
  vector<Type> jnll(3);
  jnll.setZero();

  // Put unconstrained parameters on their natural (constrained) scales
  vector<Type> kappa = exp(log_kappa);
  vector<Type> tau = exp(log_tau);
  Type sigma = exp(log_sigma);

  // Project spatial effects from mesh nodes to observation locations
  vector<Type> omega_n(N_obs);
  vector<Type> omega_w(N_obs);
  omega_n = A * spat_n;
  omega_w = A * spat_w;

  // Get density of spatial random effects, repeated for each process. Scale the
  // precision matrix by τ² to match the usual (Lindgren et al. 2011)
  // formulation.
  SparseMatrix<Type> Q_n = tau(0) * Q_spde(spde, kappa(0)) * tau(0);
  SparseMatrix<Type> Q_w = tau(1) * Q_spde(spde, kappa(1)) * tau(1);
  jnll(1) += GMRF(Q_n)(spat_n);
  jnll(2) += GMRF(Q_w)(spat_w);

  // Simulate spatial random effects using given precision matrices. Then
  // project them to the provided locations. Can't simulate new locations
  // without recomputing the A matrix, which requires the INLA package.
  SIMULATE {
    spat_n = GMRF(Q_n).simulate();
    omega_n = A * spat_n;
    spat_w = GMRF(Q_w).simulate();
    omega_w = A * spat_w;

    REPORT(spat_n);
    REPORT(spat_w);
  }

  // Get group density (n) and weight per group (w) for each observation
  vector<Type> log_n(N_obs);
  vector<Type> log_w(N_obs);
  log_n = X_n * beta_n + omega_n;
  log_w = X_w * beta_w + omega_w;

  // Convert to log-probability encounter (p), log-probability of zero catch
  // (1mp) and positive catch rate (r)
  // TODO Use a N_obs×3 matrix here instead of three vectors. Maybe transpose
  // that so that elements to be accessed together are in a column?
  vector<Type> log_p_enc(N_obs);
  vector<Type> log_p_zero(N_obs);
  vector<Type> log_r(N_obs);
  // Temporary vector to use for return of logpoislink function. See todo above
  // for probably more elegant solution.
  vector<Type> log_ppr_tmp(3);
  for (int i = 0; i < N_obs; i++) {
    // TODO VECTORIZE2_tt poislink function? Not sure if this will work if it
    // returns a vector. Probably better to pass in views to p and r here.
    log_ppr_tmp = logpoislink(log_n(i), log_w(i));
    log_p_enc(i) = log_ppr_tmp(0);
    log_p_zero(i) = log_ppr_tmp(1);
    log_r(i) = log_ppr_tmp(2);
  }

  // Calculate likelihood
  for (int i = 0; i < N_obs; i++) {
    if (catch_obs(i) == 0) {
      jnll(0) -= log_p_zero(i);
    } else {
      jnll(0) -= log_p_enc(i) +
        dlnorm(catch_obs(i), log_r(i) - sigma * sigma / Type(2.0), sigma, true);
    }
  }

  // Simulate observations. Use `p` as encounter probability, and `r` as median
  // of log Normal with logsd `sigma`.
  SIMULATE {
    vector<Type> encounter(N_obs);
    vector<Type> log_catch_median(N_obs);

    encounter = rbinom(Type(1.0), exp(log_p_enc));
    log_catch_median = log_r - sigma * sigma / Type(2.0);
    catch_obs = encounter * rlnorm(log_catch_median, sigma);

    REPORT(encounter);
    REPORT(catch_obs);
  }

  REPORT(jnll);

  return jnll.sum();
}

