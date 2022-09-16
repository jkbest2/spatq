##' Indicate which parameters should be estimated; conversely, which should be
##' mapped to their initial values. Each argument can be logical, indicating
##' whether both the numbers density, weight per group, and their
##' hyperparameters should be estimated. A list specifying whether the numbers
##' density, weight per group, and hyperparameters (potentially nested inside
##' each specific process).
##'
##' @title Specify which parameters to estimate
##' @param beta Abundance fixed effects
##' @param gamma Abundance iid random effects
##' @param omega Abundance spatial effects
##' @param epsilon Abundance spatiotemporal effects
##' @param lambda Catchability fixed effects
##' @param eta Catchability iid random effects
##' @param phi Catchability spatial effects
##' @param psi Catchability spatiotemporal effects
##' @param kappa_map Factor vector of length 8 specifying the map for kappa
##'   parameters
##' @param obs_lik Observation likelihood. Current options are \code{0}:
##'   Poisson-link log-normal and \code{1}: Tweedie
##' @return A list with elements for each random effect process and parameter in
##'   the \code{spatq} model, each with a logical indicating whether it is
##'   (\code{TRUE}) or is not (\code{FALSE}) estimated.
##' @author John Best
##' @export
specify_estimated <- function(beta = TRUE,
                              gamma = FALSE,
                              omega = FALSE,
                              epsilon = FALSE,
                              lambda = TRUE,
                              eta = FALSE,
                              phi = FALSE,
                              psi = FALSE,
                              kappa_map = NULL,
                              aniso = FALSE,
                              obs_lik = 0L) {
  estd <- list()
  estd$beta_n <- is_estd_proc("beta_n", beta, obs_lik)
  estd$beta_w <- is_estd_proc("beta_w", beta, obs_lik)
  estd$gamma_n <- is_estd_proc("gamma_n", gamma, obs_lik)
  estd$gamma_w <- is_estd_proc("gamma_w", gamma, obs_lik)
  estd$omega_n <- is_estd_proc("omega_n", omega, obs_lik)
  estd$omega_w <- is_estd_proc("omega_w", omega, obs_lik)
  estd$epsilon_n <- is_estd_proc("epsilon_n", epsilon, obs_lik)
  estd$epsilon_w <- is_estd_proc("epsilon_w", epsilon, obs_lik)

  estd$lambda_n <- is_estd_proc("lambda_n", lambda, obs_lik)
  estd$lambda_w <- is_estd_proc("lambda_w", lambda, obs_lik)
  estd$eta_n <- is_estd_proc("eta_n", eta, obs_lik)
  estd$eta_w <- is_estd_proc("eta_w", eta, obs_lik)
  estd$phi_n <- is_estd_proc("phi_n", phi, obs_lik)
  estd$phi_w <- is_estd_proc("phi_w", phi, obs_lik)
  estd$psi_n <- is_estd_proc("psi_n", psi, obs_lik)
  estd$psi_w <- is_estd_proc("psi_w", psi, obs_lik)

  spec <- list(
    gamma = gamma, omega = omega, epsilon = epsilon,
    eta = eta, phi = phi, psi = psi
  )
  estd$log_xi <- is_estd_hpar("log_xi", estd)
  estd$log_kappa <- is_estd_hpar("log_kappa", estd)
  estd$log_tau <- is_estd_hpar("log_tau", estd)
  estd$H_pars <- aniso
  estd$obs_lik <- obs_lik
  if (!is.null(kappa_map)) {
    attr(estd, "kappa_map") <- factor(kappa_map)
  }

  return(
    structure(estd,
      class = "spatq_spec"
    )
  )
}

##' Only intended to be called from \code{\link{specify_estimated}};
##' \code{proc_spec} argument should take any of the arguments of
##' \code{specify_estimated} while \code{spec} argument takes a list containing
##' all of those arguments..
##'
##' @title Check if process is estimated
##' @param parname Name of the process parameter (e.g. \code{"omega_n"}) as a string
##' @param proc_spec Logical or list at the "process" level (e.g. omega, which
##'   contains \code{omega_n} and \code{omega_w})
##' @param estd List indicating whether each random process is estimated
##' @param obs_lik Observation likelihood, current options are \code{0}:
##    Poisson-link log-normal and \code{1}: Tweedie
##' @return Logical indicating whether parameter is to be estimated
##' @author John Best
is_estd_proc <- function(parname, proc_spec, obs_lik) {
  if (obs_lik == 1 && grepl(r"(_w$)", parname)) {
    ## Don't estimate weight-per-group parameters if using a Tweedie
    ## observation likelihood
    est <- FALSE
  } else if (is.logical(proc_spec)) {
    est <- proc_spec
  } else if (is.logical(proc_spec[[parname]])) {
    est <- proc_spec[[parname]]
  } else if (is.list(proc_spec) && !(parname %in% names(proc_spec))) {
    est <- TRUE
  } else if (is.list(proc_spec[[parname]])) {
    est <- TRUE
  } else {
    ## Fallback; something went wrong
    stop("Parameters ", parname, " not specified correctly.")
  }

  return(est)
}

##' @describeIn is_estd_proc Check which hyperparameters are estimated
is_estd_hpar <- function(parname, estd) {
  if (parname == "log_xi") {
    procs <- c("gamma_n", "gamma_w", "eta_n", "eta_w")
  } else {
    procs <- rep(c("omega", "epsilon", "phi", "psi"), each = 2)
    procs <- paste0(procs, c("_n", "_w"))
  }

  unlist(estd[procs])
}
