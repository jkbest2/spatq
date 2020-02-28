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
##' @param epsilon1 Abundance spatiotemporal effects (sum-to-zero constrained)
##' @param lambda Catchability fixed effects
##' @param eta Catchability iid random effects
##' @param phi Catchability spatial effects
##' @param psi1 Catchability spatiotemporal effects (sum-to-zero constrained)
##' @return A list with elements for each random effect process and parameter in
##'   the \code{spatq} model, each with a logical indicating whether it is
##'   (\code{TRUE}) or is not (\code{FALSE}) estimated.
##' @author John Best
##' @export
specify_estimated <- function(beta = TRUE,
                              gamma = FALSE,
                              omega = FALSE,
                              epsilon1 = FALSE,
                              lambda = TRUE,
                              eta = FALSE,
                              phi = FALSE,
                              psi1 = FALSE) {
  estd <- list()
  estd$beta_n <- is_estd_proc("beta_n", beta)
  estd$beta_w <- is_estd_proc("beta_w", beta)
  estd$gamma_n <- is_estd_proc("gamma_n", gamma)
  estd$gamma_w <- is_estd_proc("gamma_w", gamma)
  estd$omega_n <- is_estd_proc("omega_n", omega)
  estd$omega_w <- is_estd_proc("omega_w", omega)
  estd$epsilon1_n <- is_estd_proc("epsilon1_n", epsilon1)
  estd$epsilon1_w <- is_estd_proc("epsilon1_w", epsilon1)

  estd$lambda_n <- is_estd_proc("lambda_n", lambda)
  estd$lambda_w <- is_estd_proc("lambda_w", lambda)
  estd$eta_n <- is_estd_proc("eta_n", eta)
  estd$eta_w <- is_estd_proc("eta_w", eta)
  estd$phi_n <- is_estd_proc("phi_n", phi)
  estd$phi_w <- is_estd_proc("phi_w", phi)
  estd$psi1_n <- is_estd_proc("psi1_n", psi1)
  estd$psi1_w <- is_estd_proc("psi1_w", psi1)

  spec <- list(gamma = gamma, omega = omega, epsilon1 = epsilon1,
               eta = eta, phi = phi, psi1 = psi1)
  estd$log_xi <- is_estd_hpar("log_xi", spec)
  estd$log_kappa <- is_estd_hpar("log_kappa", spec)
  estd$log_tau <- is_estd_hpar("log_tau", spec)

  return(estd)
}

##' Only intended to be called from \code{\link{specify_estimated}};
##' \code{proc_spec} argument should take any of the arguments of
##' \code{specify_estimated} while \code{spec} argument takes a list containing
##' all of those arguments..
##'
##' @title Check if process is estimated
##' @param par Name of the process parameter (e.g. \code{"omega_n"}) as a string
##' @param proc_spec Logical or list at the "process" level (e.g. omega, which
##'   contains \code{omega_n} and \code{omega_w})
##' @param spec List of all process-level specifications
##' @return Logical indicating whether parameter is to be estimated
##' @author John Best
is_estd_proc <- function(par, proc_spec) {
  if (is.logical(proc_spec)) {
    est <- proc_spec
  } else if (is.logical(proc_spec[[par]])) {
    est <- proc_spec[[par]]
  } else if (is.list(proc_spec) && !(par %in% names(proc_spec))) {
    est <- TRUE
  } else if (is.list(proc_spec[[par]])) {
    est <- TRUE
  } else {
    est <- NA
  }

  return(est)
}

##' @describeIn is_estd_proc Check which hyperparameters are estimated
is_estd_hpar <- function(par, spec) {
  if (par == "log_xi") {
    procs <- c("gamma", "gamma", "eta", "eta")
    pars <- c("gamma_n", "gamma_w", "eta_n", "eta_w")
  } else {
    procs <- rep(c("omega", "epsilon1", "phi", "psi1"), each = 2)
    pars <- paste0(procs, c("_n", "_w"))
  }

  ## Use NA's so a failed parse of the spec tree is obvious
  est <- rep(NA, length(pars))
  for (idx in seq_along(pars)) {
    ## Extract the process to consider (either "gamma" or "eta")
    proc_spec <- spec[[procs[idx]]]
    if (is.logical(proc_spec)) {
      ## If the process is logical, use that for hyperparameters as well
      est[idx] <- proc_spec
    } else if (is.list(proc_spec)) {
      if (par %in% names(proc_spec)) {
        ## If hyperparameters are specified for both parameters simultaneously,
        ## use those values
        est[idx] <- proc_spec[[par]]
      } else if (pars[idx] %in% names(proc_spec)) {
        ## If n and w are specified separately, extract the specific
        ## parameter of interest
        par_spec <- proc_spec[[pars[idx]]]
        if (is.logical(par_spec)) {
          ## If this value is logical, apply to log_xi as well
          est[idx] <- par_spec
        } else if (par %in% names(par_spec)) {
          ## Otherwise, take value from nested "log_xi" list element
          est[idx] <- par_spec[[par]]
        }
      }
    }
  }
  ## Return the logical vector
  return(est)
}
