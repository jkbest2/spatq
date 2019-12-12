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
  // Area swept offset for each observation
  DATA_VECTOR(area_swept);

  // Abundance fixed effects design matrix
  DATA_MATRIX(X_n);
  DATA_MATRIX(X_w);
  // Index fixed effects design matrix
  DATA_MATRIX(IX_n);
  DATA_MATRIX(IX_w);

  // Abundance random effects design matrix
  DATA_MATRIX(Z_n);
  DATA_MATRIX(Z_w);
  // Index random effects
  DATA_MATRIX(IZ_n);
  DATA_MATRIX(IZ_w);

  // Abundance projection matrices
  DATA_SPARSE_MATRIX(A_spat);      // Spatial
  DATA_SPARSE_MATRIX(A_sptemp);    // Spatiotemporal
  // Index projection matrices
  DATA_SPARSE_MATRIX(IA_spat);     // Spatial
  DATA_SPARSE_MATRIX(IA_sptemp);   // Spatiotemporal

  // Index integration weights
  DATA_VECTOR(Ih);

  // ---------------------------------------------------------------------------
  // Fixed effects design matrix
  DATA_MATRIX(R_n);
  DATA_MATRIX(R_w);

  // Random effects design matrix
  DATA_MATRIX(V_n);
  DATA_MATRIX(V_w);

  // Catchability projection matrices
  DATA_SPARSE_MATRIX(A_qspat);     // Spatial fishery-dependent
  DATA_SPARSE_MATRIX(A_qsptemp);   // Spatiotemporal fishery-dependent

  // ---------------------------------------------------------------------------
  // FEM matrices for SPDE spatial and spatiotemporal effects. Sharing the mesh
  // between effects means that this only needs to be passed once.
  DATA_STRUCT(spde, spde_t);

  // ---------------------------------------------------------------------------
  // Vector indicating which of the 12 random processes should be included.
  // Currently takes a length-6 vector, and can only switch off pairs of numbers
  // density and weight-per-group processes. Indexing is currently:
  // 0: gamma_n, gamma_w
  // 1: omega_n, omega_w
  // 2: epsilon_n, epsilon_w
  // 3: eta_n, eta_w
  // 4: phi_n, phi_w
  // 5: psi_n, psi_w
  DATA_IVECTOR(proc_switch);

  // ===========================================================================
  // PARAMETER section
  // ---------------------------------------------------------------------------
  // Abundance fixed effects
  PARAMETER_VECTOR(beta_n);        // number of fixed effects
  PARAMETER_VECTOR(beta_w);        // number of fixed effects

  // Abundance random effects
  PARAMETER_VECTOR(gamma_n);       // number of random effects
  PARAMETER_VECTOR(gamma_w);       // number of random effects

  // Abundance spatial effects
  PARAMETER_VECTOR(omega_n);       // N_vert
  PARAMETER_VECTOR(omega_w);       // N_vert

  // Abundance spatiotemporal effects
  PARAMETER_MATRIX(epsilon_n);     // N_vert × N_yrs
  PARAMETER_MATRIX(epsilon_w);     // N_vert × N_yrs

  // ---------------------------------------------------------------------------
  // Catchability fixed effects
  PARAMETER_VECTOR(lambda_n);      // number of fixed effects
  PARAMETER_VECTOR(lambda_w);      // number of fixed effects

  // Catchability random effects
  PARAMETER_VECTOR(eta_n);         // number of random effects
  PARAMETER_VECTOR(eta_w);         // number of random effects

  // Catchability spatial effects
  PARAMETER_VECTOR(phi_n);         // N_vert
  PARAMETER_VECTOR(phi_w);         // N_vert

  // Catchability spatiotemporal effects
  PARAMETER_MATRIX(psi_n);         // N_vert × N_yrs
  PARAMETER_MATRIX(psi_w);         // N_vert × N_yrs

  // ---------------------------------------------------------------------------
  // Random effects variance parameters
  PARAMETER_VECTOR(log_xi);

  // Spatial and spatiotemporal field parameters
  PARAMETER_VECTOR(log_kappa);     // 8
  PARAMETER_VECTOR(log_tau);       // 8

  // Log catch variation parameter
  PARAMETER(log_sigma);            // 1

  // ===========================================================================
  // Derived values
  // ---------------------------------------------------------------------------
  // Get number of observations
  int N_obs = catch_obs.size();
  // Get number of years
  int N_yrs = epsilon_n.cols();
  // Get number of integration locations for each index year
  int N_I = Ih.size();

  // Put unconstrained parameters on their natural (constrained) scales
  vector<Type> xi = exp(log_xi);
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
  // 8:   Abundance random effects
  // 9:   Catchability random effects
  // 10:  Observation likelihood
  vector<Type> jnll(11);
  jnll.setZero();

  // ===========================================================================
  // Abundance fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> fixef_n(N_obs);
  vector<Type> fixef_w(N_obs);
  fixef_n = X_n * beta_n;
  fixef_w = X_w * beta_w;

  // Index fixed effects
  vector<Type> Ifixef_n(N_I);
  vector<Type> Ifixef_w(N_I);
  Ifixef_n = IX_n * beta_n;
  Ifixef_w = IX_w * beta_w;

  // ===========================================================================
  // Abundance random effects
  // ---------------------------------------------------------------------------
  vector<Type> ranef_n(N_obs);
  vector<Type> ranef_w(N_obs);
  vector<Type> Iranef_n(N_I);
  vector<Type> Iranef_w(N_I);

  if (proc_switch(0)) {
    ranef_n = Z_n * gamma_n;
    ranef_w = Z_w * gamma_w;

    // Index fixed effects
    Iranef_n = IZ_n * gamma_n;
    Iranef_w = IZ_w * gamma_w;

    // Include random effects likelihood; currently only iid is implemented
    jnll(8) -= sum(dnorm(gamma_n, Type(0), xi(0), true));
    jnll(8) -= sum(dnorm(gamma_w, Type(0), xi(1), true));

    SIMULATE {
      gamma_n = rnorm(Z_n.cols(), Type(0), xi(0));
      ranef_n = Z_n * gamma_n;
      gamma_w = rnorm(Z_w.cols(), Type(0), xi(1));
      ranef_w = Z_w * gamma_w;

      REPORT(gamma_n);
      REPORT(gamma_w);
    }
  } else {
    ranef_n.setZero();
    ranef_w.setZero();
    Iranef_n.setZero();
    Iranef_w.setZero();
  }

  // ===========================================================================
  // Abundance spatial effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> spat_n(N_obs);
  vector<Type> spat_w(N_obs);
  vector<Type> Ispat_n(N_I);
  vector<Type> Ispat_w(N_I);

  if (proc_switch(1)) {
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

    // Index spatial effects
    Ispat_n = IA_spat * omega_n;
    Ispat_w = IA_spat * omega_w;
  } else {
    spat_n.setZero();
    spat_w.setZero();
    Ispat_n.setZero();
    Ispat_w.setZero();
  }

  // ===========================================================================
  // Abundance spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  // TODO Implement (optional?) sum-to-zero constraint
  vector<Type> sptemp_n(N_obs);
  vector<Type> sptemp_w(N_obs);
  vector<Type> Isptemp_n(N_I);
  vector<Type> Isptemp_w(N_I);

  if (proc_switch(2)) {
    sptemp_n = A_sptemp * epsilon_n.value();
    sptemp_w = A_sptemp * epsilon_w.value();

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

      REPORT(epsilon_n);
      REPORT(epsilon_w);
    }

    // Index spatiotemporal effects
    Isptemp_n = IA_sptemp * epsilon_n.value();
    Isptemp_w = IA_sptemp * epsilon_w.value();
  } else {
    sptemp_n.setZero();
    sptemp_w.setZero();
    Isptemp_n.setZero();
    Isptemp_w.setZero();
  }

  // ===========================================================================
  // Catchability fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> qfixef_n(N_obs);
  vector<Type> qfixef_w(N_obs);
  qfixef_n = R_n * lambda_n;
  qfixef_w = R_w * lambda_w;

  // ===========================================================================
  // Catchability random effects
  // ---------------------------------------------------------------------------
  vector<Type> qranef_n(N_obs);
  vector<Type> qranef_w(N_obs);

  if (proc_switch(3)) {
    qranef_n = V_n * eta_n;
    qranef_w = V_w * eta_w;

    // Include random effects likelihood; currently only iid is implemented
    jnll(9) -= sum(dnorm(eta_n, Type(0), xi(2), true));
    jnll(9) -= sum(dnorm(eta_w, Type(0), xi(3), true));

    SIMULATE {
      eta_n = rnorm(V_n.cols(), Type(0), xi(2));
      qranef_n = V_n * eta_n;
      eta_w = rnorm(V_w.cols(), Type(0), xi(3));
      qranef_w = V_w * eta_w;

      REPORT(eta_n);
      REPORT(eta_w);
    }
  } else {
    qranef_n.setZero();
    qranef_w.setZero();
  }

  // ===========================================================================
  // Catchability spatial effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> qspat_n(N_obs);
  vector<Type> qspat_w(N_obs);

  if (proc_switch(4)) {
    qspat_n = A_qspat * phi_n;
    qspat_w = A_qspat * phi_w;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_ph = tau(4) * Q_spde(spde, kappa(4)) * tau(4);
    SparseMatrix<Type> Q_w_ph = tau(5) * Q_spde(spde, kappa(5)) * tau(5);
    jnll(4) += GMRF(Q_n_ph)(phi_n);
    jnll(5) += GMRF(Q_w_ph)(phi_w);

    // Simulate spatial random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      phi_n = GMRF(Q_n_ph).simulate();
      qspat_n = A_qspat * phi_n;
      phi_w = GMRF(Q_w_ph).simulate();
      qspat_w = A_qspat * phi_w;

      REPORT(phi_n);
      REPORT(phi_w);
    }
  } else {
    qspat_n.setZero();
    qspat_w.setZero();
  }

  // ===========================================================================
  // Catchability spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  // TODO Implement (optional?) sum-to-zero constraint
  vector<Type> qsptemp_n(N_obs);
  vector<Type> qsptemp_w(N_obs);

  if (proc_switch(5)) {
    qsptemp_n = A_qsptemp * psi_n.value();
    qsptemp_w = A_qsptemp * psi_w.value();

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_ps = tau(6) * Q_spde(spde, kappa(6)) * tau(6);
    GMRF_t<Type> gmrf_n_ps(Q_n_ps);
    SparseMatrix<Type> Q_w_ps = tau(7) * Q_spde(spde, kappa(7)) * tau(7);
    GMRF_t<Type> gmrf_w_ps(Q_w_ps);
    for (int yr = 0; yr < N_yrs; yr++) {
      jnll(6) += gmrf_n_ps(psi_n.col(yr));
      jnll(7) += gmrf_w_ps(psi_w.col(yr));
    }

    // Simulate spatiotemporal random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      for (int yr = 0; yr < N_yrs; yr++) {
        psi_n.col(yr) = gmrf_n_ps.simulate();
        psi_w.col(yr) = gmrf_w_ps.simulate();
      }
      qsptemp_n = A_qsptemp * psi_n.value();
      qsptemp_w = A_qsptemp * psi_w.value();

      REPORT(psi_n);
      REPORT(psi_w);
    }
  } else {
    qsptemp_n.setZero();
    qsptemp_w.setZero();
  }

  // ===========================================================================
  // Calculate linear predictor
  // ---------------------------------------------------------------------------
  // Get group density (n) and weight per group (w) for each observation
  vector<Type> log_n(N_obs);
  vector<Type> log_w(N_obs);
  log_n = fixef_n + ranef_n + spat_n + sptemp_n +
    qfixef_n + qranef_n + qspat_n + qsptemp_n;
  log_w = fixef_w + ranef_w + spat_w + sptemp_w +
    qfixef_w + qranef_w + qspat_w + qsptemp_w;

  // Index linear predictor
  vector<Type> Ilog_n(N_I);
  vector<Type> Ilog_w(N_I);
  Ilog_n = Ifixef_n + Iranef_n + Ispat_n + Isptemp_n;
  Ilog_w = Ifixef_w + Iranef_w + Ispat_w + Isptemp_w;

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
    log_ppr_tmp = logpoislink(log_n(i), log_w(i), area_swept(i));
    log_p_enc(i) = log_ppr_tmp(0);
    log_p_zero(i) = log_ppr_tmp(1);
    log_r(i) = log_ppr_tmp(2);
  }

  // ===========================================================================
  // Observation likelihood
  // ---------------------------------------------------------------------------
  for (int i = 0; i < N_obs; i++) {
    if (catch_obs(i) == 0) {
      jnll(10) -= log_p_zero(i);
    } else {
      jnll(10) -= log_p_enc(i) +
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
  // Index of abundance
  // ---------------------------------------------------------------------------
  vector<Type> Index(N_yrs);
  Index.setZero();
  int i0 = 0;
  int i1 = 0;

  int N_I_per_yr = N_I / N_yrs;

  for (int yr = 0; yr < N_yrs; yr++) {
    i0 = yr * N_I_per_yr;
    i1 = (yr + 1) * N_I_per_yr;
    for (int i = i0; i < i1; i++) {
      // Should prevent some numerical issues; cribbed from VAST index
      // calculation.
      Index(yr) += Ih(i) * exp(Ilog_n(i)) * exp(Ilog_w(i));
    }
  }

  // ===========================================================================
  // Reports
  // ---------------------------------------------------------------------------
  vector<Type> rho_sp;
  vector<Type> sigma_sp;
  rho_sp = sqrt(8) / kappa;
  sigma_sp = 1 / (kappa * tau * 2 * sqrt(PI));

  REPORT(jnll);
  REPORT(Ilog_n);
  REPORT(Ilog_w);
  REPORT(rho_sp);
  REPORT(sigma_sp);

  ADREPORT(Index);

  return jnll.sum();
}

