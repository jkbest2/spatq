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
  // Vector indicating which of the 12 random processes should be included in
  // the log-likelihood calculation. Takes a length-12 vector, and can switch
  // off individual random processes. Indexing is currently:
  // 0: gamma_n
  // 1: gamma_w
  // 2: omega_n
  // 3: omega_w
  // 4: epsilon_n
  // 5: epsilon_w
  // 6: eta_n
  // 7: eta_w
  // 8: phi_n
  // 9: phi_w
  // 10: psi_n
  // 11: psi_w
  DATA_IVECTOR(proc_switch);

  // ---------------------------------------------------------------------------
  // Flags to control GMRF normalization and return early for normalization
  DATA_INTEGER(norm_flag);
  DATA_INTEGER(incl_data);

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
  PARAMETER_VECTOR(epsilon_n);     // (N_vert × N_yrs)
  PARAMETER_VECTOR(epsilon_w);     // (N_vert × N_yrs)

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
  PARAMETER_VECTOR(psi_n);         // (N_vert × N_yrs)
  PARAMETER_VECTOR(psi_w);         // (N_vert × N_yrs)

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
  // Get number of vertices in mesh
  int N_vert = omega_n.size();
  // Get number of observations
  int N_obs = catch_obs.size();
  // Get number of years
  int N_yrs = epsilon_n.size() / N_vert;
  // Get number of integration locations for each index year
  int N_I = Ih.size();
  // Convert norm_flag incl_data to boolean
  bool nrmlz = bool(norm_flag);
  bool incl_datalik = bool(incl_data);

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

  // gamma_n
  if (proc_switch(0)) {
    ranef_n = Z_n * gamma_n;
    // Index fixed effects
    Iranef_n = IZ_n * gamma_n;

    // Include random effects likelihood; currently only iid is implemented
    jnll(8) -= sum(dnorm(gamma_n, Type(0), xi(0), true));

    SIMULATE {
      gamma_n = rnorm(Z_n.cols(), Type(0), xi(0));
      ranef_n = Z_n * gamma_n;

      REPORT(gamma_n);
    }
  } else {
    ranef_n.setZero();
    Iranef_n.setZero();
  }

  // gamma_w
  if (proc_switch(1)) {
    ranef_w = Z_w * gamma_w;

    // Index fixed effects
    Iranef_w = IZ_w * gamma_w;

    // Include random effects likelihood; currently only iid is implemented
    jnll(8) -= sum(dnorm(gamma_w, Type(0), xi(1), true));

    SIMULATE {
      gamma_w = rnorm(Z_w.cols(), Type(0), xi(1));
      ranef_w = Z_w * gamma_w;

      REPORT(gamma_w);
    }
  } else {
    ranef_w.setZero();
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

  // omega_n
  if (proc_switch(2)) {
    spat_n = A_spat * omega_n;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_om = Q_spde(spde, kappa(0));
    SCALE_t<GMRF_t<Type>> gmrf_n_om = SCALE(GMRF(Q_n_om, nrmlz), tau(0));
    jnll(0) += gmrf_n_om(omega_n);

    // Simulate spatial random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      gmrf_n_om.simulate(omega_n);
      spat_n = A_spat * omega_n;

      REPORT(omega_n);
    }

    // Index spatial effects
    Ispat_n = IA_spat * omega_n;
  } else {
    spat_n.setZero();
    Ispat_n.setZero();
  }

  // omega_w
  if (proc_switch(3)) {
    spat_w = A_spat * omega_w;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_w_om = Q_spde(spde, kappa(1));
    SCALE_t<GMRF_t<Type>> gmrf_w_om = SCALE(GMRF(Q_w_om, nrmlz), tau(1));
    jnll(1) += gmrf_w_om(omega_w);

    // Simulate spatial random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      gmrf_w_om.simulate(omega_w);
      spat_w = A_spat * omega_w;

      REPORT(omega_w);
    }

    // Index spatial effects
    Ispat_w = IA_spat * omega_w;
  } else {
    spat_w.setZero();
    Ispat_w.setZero();
  }

  // ===========================================================================
  // Abundance spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> sptemp_n(N_obs);
  vector<Type> sptemp_w(N_obs);
  vector<Type> Isptemp_n(N_I);
  vector<Type> Isptemp_w(N_I);

  // epsilon_n
  if (proc_switch(4)) {
    sptemp_n = A_sptemp * epsilon_n;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_ep = Q_spde(spde, kappa(2));
    SCALE_t<GMRF_t<Type>> gmrf_n_ep = SCALE(GMRF(Q_n_ep, nrmlz), tau(2));

    for (int yr = 0; yr < N_yrs; yr++) {
      int i0 = N_vert * yr;
      jnll(2) += gmrf_n_ep(epsilon_n.segment(i0, N_vert));
    }

    // Simulate spatiotemporal random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      // Declare temporary vector to hold simulated effects. Prevents "non-const
      // lvalue" error.
      vector<Type> sim_temp_n(N_vert);

      for (int yr = 0; yr < N_yrs; yr++) {
        int i0 = N_vert * yr;
        gmrf_n_ep.simulate(sim_temp_n);
        epsilon_n.segment(i0, N_vert) = sim_temp_n;
      }

      sptemp_n = A_sptemp * epsilon_n;

      REPORT(epsilon_n);
    }

    // Index spatiotemporal effects
    Isptemp_n = IA_sptemp * epsilon_n;
  } else {
    epsilon_n.setZero();
    sptemp_n.setZero();
    Isptemp_n.setZero();
  }

  // epsilon_w
  if (proc_switch(5)) {
    sptemp_w = A_sptemp * epsilon_w;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_w_ep = Q_spde(spde, kappa(3));
    SCALE_t<GMRF_t<Type>> gmrf_w_ep = SCALE(GMRF(Q_w_ep, nrmlz), tau(3));

    for (int yr = 0; yr < N_yrs; yr++) {
      int i0 = N_vert * yr;
      jnll(3) += gmrf_w_ep(epsilon_w.segment(i0, N_vert));
    }

    // Simulate spatiotemporal random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      // Declare temporary vector to hold simulated effects. Prevents "non-const
      // lvalue" error.
      vector<Type> sim_temp_w(N_vert);

      for (int yr = 0; yr < N_yrs; yr++) {
        int i0 = N_vert * yr;
        gmrf_w_ep.simulate(sim_temp_w);
        epsilon_w.segment(i0, N_vert) = sim_temp_w;
      }

      sptemp_w = A_sptemp * epsilon_w;

      REPORT(epsilon_w);
    }

    // Index spatiotemporal effects
    Isptemp_w = IA_sptemp * epsilon_w;
  } else {
    epsilon_w.setZero();
    sptemp_w.setZero();
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

  // eta_n
  if (proc_switch(6)) {
    qranef_n = V_n * eta_n;

    // Include random effects likelihood; currently only iid is implemented
    jnll(9) -= sum(dnorm(eta_n, Type(0), xi(2), true));

    SIMULATE {
      eta_n = rnorm(V_n.cols(), Type(0), xi(2));
      qranef_n = V_n * eta_n;

      REPORT(eta_n);
    }
  } else {
    qranef_n.setZero();
  }

  // eta_w
  if (proc_switch(7)) {
    qranef_w = V_w * eta_w;

    // Include random effects likelihood; currently only iid is implemented
    jnll(9) -= sum(dnorm(eta_w, Type(0), xi(3), true));

    SIMULATE {
      eta_w = rnorm(V_w.cols(), Type(0), xi(3));
      qranef_w = V_w * eta_w;

      REPORT(eta_w);
    }
  } else {
    qranef_w.setZero();
  }

  // ===========================================================================
  // Catchability spatial effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> qspat_n(N_obs);
  vector<Type> qspat_w(N_obs);

  // phi_n
  if (proc_switch(8)) {
    qspat_n = A_qspat * phi_n;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_ph = Q_spde(spde, kappa(4));
    SCALE_t<GMRF_t<Type>> gmrf_n_ph = SCALE(GMRF(Q_n_ph, nrmlz), tau(4));
    jnll(4) += gmrf_n_ph(phi_n);

    // Simulate spatial random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      gmrf_n_ph.simulate(phi_n);
      qspat_n = A_qspat * phi_n;

      REPORT(phi_n);
    }
  } else {
    qspat_n.setZero();
  }

  // phi_w
  if (proc_switch(9)) {
    qspat_w = A_qspat * phi_w;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_w_ph = Q_spde(spde, kappa(5));
    SCALE_t<GMRF_t<Type>> gmrf_w_ph = SCALE(GMRF(Q_w_ph, nrmlz), tau(5));
    jnll(5) += gmrf_w_ph(phi_w);

    // Simulate spatial random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      gmrf_w_ph.simulate(phi_w);
      qspat_w = A_qspat * phi_w;

      REPORT(phi_w);
    }
  } else {
    qspat_w.setZero();
  }

  // ===========================================================================
  // Catchability spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> qsptemp_n(N_obs);
  vector<Type> qsptemp_w(N_obs);

  // psi_n
  if (proc_switch(10)) {
    qsptemp_n = A_qsptemp * psi_n;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_n_ps = Q_spde(spde, kappa(6));
    SCALE_t<GMRF_t<Type>> gmrf_n_ps = SCALE(GMRF(Q_n_ps, nrmlz), tau(6));

    for (int yr = 0; yr < N_yrs; yr++) {
      int i0 = N_vert * yr;
      jnll(6) += gmrf_n_ps(psi_n.segment(i0, N_vert));
    }

    // Simulate spatiotemporal random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      // Declare temporary vector to hold simulated effects. Prevents "non-const
      // lvalue" error.
      vector<Type> sim_temp_n(N_vert);

      for (int yr = 0; yr < N_yrs; yr++) {
        int i0 = N_vert * yr;
        gmrf_n_ps.simulate(sim_temp_n);
        psi_n.segment(i0, N_vert) = sim_temp_n;
      }

      qsptemp_n = A_qsptemp * psi_n;

      REPORT(psi_n);
    }
  } else {
    qsptemp_n.setZero();
  }

  // psi_w
  if (proc_switch(11)) {
    qsptemp_w = A_qsptemp * psi_w;

    // Get density of spatial random effects, repeated for each process. Scale the
    // precision matrix by τ² to match the usual (Lindgren et al. 2011)
    // formulation.
    SparseMatrix<Type> Q_w_ps = Q_spde(spde, kappa(7));
    SCALE_t<GMRF_t<Type>> gmrf_w_ps = SCALE(GMRF(Q_w_ps, nrmlz), tau(7));

    for (int yr = 0; yr < N_yrs; yr++) {
      int i0 = N_vert * yr;
      jnll(7) += gmrf_w_ps(psi_w.segment(i0, N_vert));
    }

    // Simulate spatiotemporal random effects using given precision matrices. Then
    // project them to the provided locations. Can't simulate new locations
    // without recomputing the A matrix, which requires the INLA package.
    SIMULATE {
      // Declare temporary vector to hold simulated effects. Prevents "non-const
      // lvalue" error.
      vector<Type> sim_temp_w(N_vert);

      for (int yr = 0; yr < N_yrs; yr++) {
        int i0 = N_vert * yr;
        gmrf_w_ps.simulate(sim_temp_w);
        psi_w.segment(i0, N_vert) = sim_temp_w;
      }

      qsptemp_w = A_qsptemp * psi_w;

      REPORT(psi_w);
    }
  } else {
    qsptemp_w.setZero();
  }

  // If normalizing externally, return the sum of the negative log-liklihood
  // early (the data likelihood and any map'd random effects are set to zero in
  // jnll, so the sum of the vector works).
  if (!incl_datalik) {
    return jnll.sum();
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
  sigma_sp = tau / (kappa * 2 * sqrt(PI));

  REPORT(jnll);
  REPORT(Ilog_n);
  REPORT(Ilog_w);
  REPORT(rho_sp);
  REPORT(sigma_sp);

  ADREPORT(Index);

  return jnll.sum();
}

