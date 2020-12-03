##' Estimates numbers density fixed effects through a binomial regression with
##' the complementary log-log link function. Log numbers density and log
##' encounter probability provide an offset for a log normal regression on
##' positive catches, hopefully providing reasonable starting parameter values
##' for weight per group covariates. Log numbers density is NOT currently
##' bias-corrected, but should still work well as starting values.
##'
##' @title Initial estimates of fixed effect parameters
##' @param data A data list, as produced by \code{prepare_data}
##' @return A list with estimates of \code{beta_n}, \code{beta_w},
##'   \code{lambda_n}, and \code{lambda_w}
##' @author John Best
##' @export
init_fixef <- function(data) {
  ## Don't include catchability effects if they aren't used in the model (e.g.
  ## if estimating index using single survey vessel).
  q_used <- any(data$R_n != 0)
  if (q_used) {
    df_n <- tibble::tibble(enc = data$catch_obs > 0,
                           X_n = data$X_n,
                           R_n = data$R_n)
  } else {
    df_n <- tibble::tibble(enc = data$catch_obs > 0,
                           X_n = data$X_n)
  }
  ## Estimate each group density fixed effect with a simple GLM, using the
  ## complementary log-log link for encounter. This corresponds directly to
  ## group density in the Poisson link model.
  mod_n <- stats::glm(enc ~ 0 + ., data = df_n,
                      family = stats::binomial(cloglog))
  est_n <- stats::coef(mod_n)

  ## Using the estimated group density and probability of encounter, calculate
  ## the log mean weight per group, conditional on a positive catch.
  enc <- df_n$enc
  if (q_used) {
    df_w <- tibble::tibble(catch_obs = data$catch_obs[enc],
                           X_w = data$X_w[enc, ],
                           R_w = data$R_w[enc, ])
  } else {
    df_w <- tibble::tibble(catch_obs = data$catch_obs[enc],
                           X_w = data$X_w[enc, ])
  }
  ## Log-positive catch rate is log(n) - log(p) + log(w), so the first two are
  ## used as an offset here.
  offset_w <- stats::predict(mod_n) -
    log(stats::predict(mod_n, type = "response"))
  offset_w <- offset_w[data$catch_obs > 0]
  ## mod_w <- glm(catch_obs ~ 0 + ., data = df_w, offset = offset_w,
  ##              family = gaussian(log))
  mod_w <- stats::lm(log(catch_obs) ~ 0 + ., data = df_w, offset = offset_w)
  ## Not going to worry about the bias correction for these initial values
  est_w <- stats::coef(mod_w)

  init <- list(beta_n = est_n[seq_len(ncol(data$X_n))],
               beta_w = est_w[seq_len(ncol(data$X_w))],
               lambda_n = ifelse(q_used, est_n[-seq_len(ncol(data$X_n))], 0),
               lambda_w = ifelse(q_used, est_w[-seq_len(ncol(data$X_w))], 0))
  lapply(init, unname)
}

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
                     c(log(1), 1.5)) # log-dispersion and shift-logit shape for Tweedie
  pars <- list(beta_n = pars_data(data$X_n),
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
               obs_lik_pars = init_olp)
  if (init_fixef) {
    init_est <- init_fixef(data)
    pars$beta_n <- init_est$beta_n
    pars$beta_w <- init_est$beta_w
    pars$lambda_n <- init_est$lambda_n
    pars$lambda_w <- init_est$lambda_w
  }
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
  spat_pars <- c("omega_n", "omega_w",
                 "epsilon_n", "epsilon_w",
                 "phi_n", "phi_w",
                 "psi_n", "psi_w")
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
  re_pars <- c("gamma_n", "gamma_w",
               "eta_n", "eta_w")
  which(par_name == re_pars)
}
