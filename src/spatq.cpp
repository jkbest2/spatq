#define TMB_LIB_INIT R_init_spatq
#include <TMB.hpp>
#include "../inst/include/spatq.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // ===========================================================================
  // DATA section
  // ---------------------------------------------------------------------------
  // Vector of observed catches; zero or positive
  DATA_VECTOR(catch_obs);

  // Abundnce fixed effects design matrix
  DATA_MATRIX(X_n);
  DATA_MATRIX(X_w);

  // Abundance projection matrices
  DATA_SPARSE_MATRIX(A_spat);      // Spatial
  DATA_SPARSE_MATRIX(A_sptemp);    // Spatiotemporal

  // ---------------------------------------------------------------------------
  // FEM matrices for SPDE spatial and spatiotemporal effects. Sharing the mesh
  // between effects means that this only needs to be passed once.
  DATA_STRUCT(spde, spde_t);

  // Need up to four projection matrices; one for spatial effects, one for
  // spatiotemporal effects, and the same two but that only apply the effects to
  // fishery-dependent observations

  // ===========================================================================
  // PARAMETER section
  // ---------------------------------------------------------------------------
  // Abundance fixed effects
  PARAMETER_VECTOR(beta_n);       // number of fixed effects
  PARAMETER_VECTOR(beta_w);       // number of fixed effects

  // Abundance spatial effects
  PARAMETER_VECTOR(omega_n);      // N_vert
  PARAMETER_VECTOR(omega_w);      // N_vert

  // Abundance spatiotemporal effects
  PARAMETER_MATRIX(epsilon_n);    // N_vert × N_yrs
  PARAMETER_MATRIX(epsilon_w);    // N_vert × N_yrs

  // ---------------------------------------------------------------------------
  // Spatial and spatiotemporal field parameters
  PARAMETER_VECTOR(log_kappa);    // 8
  PARAMETER_VECTOR(log_tau);      // 8

  // Log catch variation parameter
  PARAMETER(log_sigma);           // 1

  // ===========================================================================
  // Derived values
  // ---------------------------------------------------------------------------
  // Get number of observations
  int N_obs = catch_obs.size();
  // Get number of years
  int N_yrs = epsilon_n.cols();

  // Put unconstrained parameters on their natural (constrained) scales
  vector<Type> kappa = exp(log_kappa);
  vector<Type> tau = exp(log_tau);
  Type sigma = exp(log_sigma);

  // ===========================================================================
  // Log-likelihood accumulator
  // ---------------------------------------------------------------------------
  // Initialize joint negative log-likelihood accumulator; accumulate spatial
  // and spatiotemporal effects in indices corresponding to indices of kappa and
  // tau.
  // 0,1: Abundance spatial
  // 2,3: Abundance spatiotemporal
  // 4,5: Catchability spatial
  // 6,7: Catchability spatiotemporal
  // 8: Observation likelihood
  vector<Type> jnll(9);
  jnll.setZero();

  // ===========================================================================
  // Abundance fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> fixef_n(N_obs);
  vector<Type> fixef_w(N_obs);
  fixef_n = X_n * beta_n;
  fixef_w = X_w * beta_w;

  // ===========================================================================
  // Abundance spatial effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> spat_n(N_obs);
  vector<Type> spat_w(N_obs);
  spat_n = A_spat * omega_n;
  spat_w = A_spat * omega_w;

  // Get density of spatial random effects, repeated for each process. Scale the
  // precision matrix by τ² to match the usual (Lindgren et al. 2011)
  // formulation.
  SparseMatrix<Type> Q_n_om = tau(0) * Q_spde(spde, kappa(0)) * tau(0);
  SparseMatrix<Type> Q_w_om = tau(1) * Q_spde(spde, kappa(1)) * tau(1);
  jnll(0) += GMRF(Q_n_om)(omega_n);
  jnll(1) += GMRF(Q_w_om)(omega_w);

  // Simulate spatial random effects using given precision matrices. Then
  // project them to the provided locations. Can't simulate new locations
  // without recomputing the A matrix, which requires the INLA package.
  SIMULATE {
    omega_n = GMRF(Q_n_om).simulate();
    spat_n = A_spat * omega_n;
    omega_w = GMRF(Q_w_om).simulate();
    spat_w = A_spat * omega_w;

    REPORT(omega_n);
    REPORT(omega_w);
  }

  // ===========================================================================
  // Abundance spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  // TODO Implement (optional?) sum-to-zero constraint
  vector<Type> sptemp_n(N_obs);
  vector<Type> sptemp_w(N_obs);
  sptemp_n = A_sptemp * epsilon_n;
  sptemp_w = A_sptemp * epsilon_w;

  // Get density of spatial random effects, repeated for each process. Scale the
  // precision matrix by τ² to match the usual (Lindgren et al. 2011)
  // formulation.
  SparseMatrix<Type> Q_n_ep = tau(2) * Q_spde(spde, kappa(2)) * tau(2);
  GMRF_t<Type> gmrf_n_ep(Q_n_ep);
  SparseMatrix<Type> Q_w_ep = tau(3) * Q_spde(spde, kappa(3)) * tau(3);
  GMRF_t<Type> gmrf_w_ep(Q_w_ep);
  for (int yr = 0; yr < N_yrs; yr++) {
    jnll(2) += gmrf_n_ep(epsilon_n.col(yr));
    jnll(3) += gmrf_w_ep(epsilon_w.col(yr));
  }

  // Simulate spatiotemporal random effects using given precision matrices. Then
  // project them to the provided locations. Can't simulate new locations
  // without recomputing the A matrix, which requires the INLA package.
  SIMULATE {
    for (int yr = 0; yr < N_yrs; yr++) {
      // TODO Figure out if these can be `gmrf_n_ep.simulate(epsilon_n)`
      epsilon_n.col(yr) = gmrf_n_ep.simulate();
      epsilon_w.col(yr) = gmrf_w_ep.simulate();
    }
    sptemp_n = A_sptemp * epsilon_n.value();
    sptemp_w = A_sptemp * epsilon_w.value();

    REPORT(omega_n);
    REPORT(omega_w);
  }

  // ===========================================================================
  // Calculate linear predictor
  // ---------------------------------------------------------------------------
  // Get group density (n) and weight per group (w) for each observation
  vector<Type> log_n(N_obs);
  vector<Type> log_w(N_obs);
  log_n = fixef_n + spat_n + sptemp_n;
  log_w = fixef_w + spat_w + sptemp_w;

  // ===========================================================================
  // Apply link function
  // ---------------------------------------------------------------------------
  // Get group density (n) and weight per group (w) for each observation
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

  // ===========================================================================
  // Observation likelihood
  // ---------------------------------------------------------------------------
  for (int i = 0; i < N_obs; i++) {
    if (catch_obs(i) == 0) {
      jnll(8) -= log_p_zero(i);
    } else {
      jnll(8) -= log_p_enc(i) +
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

  // ===========================================================================
  // Reports
  // ---------------------------------------------------------------------------
  REPORT(jnll);

  return jnll.sum();
}

