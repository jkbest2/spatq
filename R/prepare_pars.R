##' Prepare a list of parameters of appropriate dimension. Otherwise
##' uninitialized parameter starting values are set to zeros except kappa and
##' tau, which use a correlation range of 50 (half the domain) and marginal
##' variance of 1, tranformed to the kappa and tau parameterization. Initial
##' estimates of fixed effects can be estimated using GLMs with the
##' \code{init_fixef = TRUE} argument. This is recommended.
##'
##' @title Prepare parameter list
##' @param data Data list, as from \code{prepare_data}
##' @param mesh SPDE mesh
##' @param init_fixef Use GLMs to estimate fixed effect parameters for abundance
##'   and catchability in both numbers density and weight per group processes
##' @return List with scalar, vector, or matrices of zeros as starting values
##'   for optimization
##' @author John Best
##' @export
prepare_pars <- function(data, mesh, init_fixef = FALSE) {
  T <- attr(data, "T")
  init_olp <- switch(data$obs_lik + 1,
    log(1), # log-dispersion for Poisson-link log-normal
    c(log(1), 1.5) # dispersion and shape parameters for Tweedie
  ) # log-dispersion and shift-logit shape for Tweedie
  pars <- list(
    beta_n = pars_data(data$X_n),
    beta_w = pars_data(data$X_w),
    gamma_n = pars_data(data$Z_n),
    gamma_w = pars_data(data$Z_w),
    omega_n = pars_data(mesh),
    omega_w = pars_data(mesh),
    epsilon_n = pars_data(mesh, T),
    epsilon_w = pars_data(mesh, T),
    lambda_n = pars_data(data$R_n),
    lambda_w = pars_data(data$R_w),
    eta_n = pars_data(data$V_n),
    eta_w = pars_data(data$V_w),
    phi_n = pars_data(mesh),
    phi_w = pars_data(mesh),
    psi_n = pars_data(mesh, T),
    psi_w = pars_data(mesh, T),
    log_xi = rep(0.0, 4L),
    log_kappa = rep(log(pars_kappa(50)), 8),
    log_tau = rep(log(pars_tau(1.0, 50)), 8),
    H_pars = c(log(1.0), 0.0),
    obs_lik_pars = init_olp
  )

  ## Add attribute noting whether `lambda_n` and `lambda_w` are map'd or not.
  ## Useful for consistency checks later.
  attr(pars, "map_lambda") <- all(data$R_n == 0)
  pars
}

##' Spatial/spatiotemporal parameters (kappa and tau) are stored in a length-8
##' vector. This function takes a name and returns the corresponding index.
##'
##' @title Spatial parameter index
##' @param par_name Name of the spatial/spatiotemporal parameter (e.g.
##'   \code{omega_n}, \code{phi_w} &c.)
##' @return The index of the effect in the parameter vectors
##' @author John Best
##' @export
spat_par_idx <- function(par_name) {
  spat_pars <- c(
    "omega_n", "omega_w",
    "epsilon_n", "epsilon_w",
    "phi_n", "phi_w",
    "psi_n", "psi_w"
  )
  which(par_name == spat_pars)
}

##' The variance of general random effects are stored in a length-4 vector. This
##' function takes a name and returns the corresponding index.
##'
##' @title General random effects parameter index
##' @param par_name Name of the spatial/spatiotemporal parameter (e.g.
##'   \code{gamm_n}, \code{eta_w} &c.)
##' @return The index of the effect in the parameter vectors
##' @author John Best
##' @export
re_par_idx <- function(par_name) {
  re_pars <- c(
    "gamma_n", "gamma_w",
    "eta_n", "eta_w"
  )
  which(par_name == re_pars)
}
