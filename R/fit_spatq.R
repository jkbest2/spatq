##' Uses \code{\link[stats]{optim}} to find an optimum. If a previous fit is not
##' provided, the first fit uses the conjugate gradient method to improve
##' quickly. After that, the BFGS method is used to refine the optimization.
##' Termination tolerances and control parameters for \code{optim} can be
##' provided through \code{\link{spatq_optcontrol}}.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param fit Previous fit to use as starting values
##' @param optcontrol a
##' @param bias.correct Use bias correction for \code{\link[TMB]{sdreport}}?
##' @param getJointPrecision Return the joint fixed and random effect precision
##'   matrix from \code{\link[TMB]{sdreport}}?
##' @param ... additional arguments to pass to \code{\link[TMB]{sdreport}}
##' @return An optimization object, report list, or sdreport list
##' @author John Best
##' @importFrom stats optim
##' @export
fit_spatq <- function(obj, fit = NULL, optcontrol = spatq_optcontrol()) {
  if (is.null(fit)) {
    fit <- init_spatq_fit(obj)
  }
  while (!terminate_opt(fit, optcontrol)) {
    fit1 <- optim(fit$par, obj$fn, obj$gr, method = "BFGS",
                  control = optcontrol$bfgs_control)
    fit <- attach_optdiags(fit1, fit, obj)
  }
  fit
}

##' Attach optimization diagnostics comparing a previous optimization result to
##' a new optimization result.
##'
##' @title Attach optimization diagnostics
##' @param newfit new optimization result
##' @param oldfit previous optimization result
##' @param obj TMB objective function
##' @return A \code{spatq_fit} object
##' @author John K Best
##' @export
attach_optdiags <- function(newfit, oldfit = NULL, obj) {
  if (is.null(oldfit)) {
    oldfit <- init_spatq_fit(obj)
  }
  newfit$grad <- obj$gr(newfit$par)
  newfit$dgrad <- oldfit$grad - newfit$grad
  newfit$dparrel <- (oldfit$par - newfit$par) / oldfit$par
  newfit$dobjrel <- (oldfit$value - newfit$value) / oldfit$value
  newfit$nfit <- oldfit$nfit + 1
  newfit$totcounts <- oldfit$counts + newfit$counts
  class(newfit) <- "spatq_fit"
  newfit
}

##' @describeIn attach_optdiags Initialize a spatq_fit object
##' @export
init_spatq_fit <- function(obj) {
  ## Fill in with values that make sense or else values that won't cause
  ## optimization to terminate prematurely
  structure(
    list(value = obj$fn(obj$par),
         par = obj$par,
         grad = obj$gr(obj$par),
         dparrel = rep_len(Inf, length(obj$par)),
         dobjrel = Inf,
         nfit = 0,
         totcounts = c("function" = 0, "gradient" = 0)),
    class = "spatq_fit")
}

##' Sequentially tests whether maximum number of optimizations has been reached,
##' if the maximum gradient component is less than \code{optcontrol$grtol},
##' relative parameter value changes between optimizations is less than
##' \code{optcontrol$dparrtol}, or relative change in the objective function
##' value is less than \code{optcontrol$dobjrtol}.
##'
##' @title Test whether optimization should be stopped
##' @param fit \code{spatq_fit} object from \code{\link{attach_optdiags}}
##' @param optcontrol \code{\link{spatq_optcontrol}} object
##' @return boolean indicating termination (TRUE) or not (FALSE)
##' @author John K Best
##' @export
terminate_opt <- function(fit = NULL, optcontrol) {
  term <- FALSE
  if (is.null(fit)) {
    term <- FALSE
  } else if (fit$nfit >= optcontrol$maxopts) {
    term <- TRUE
  ## } else if (max(abs(fit$grad)) < optcontrol$grtol) {
  ##   term <- TRUE
  ## } else if (max(abs(fit$dparrel)) < optcontrol$dparrtol) {
  ##   term <- TRUE
  ## } else if (fit$dobjrel < optcontrol$dobjrtol) {
  ##   term <- TRUE
  }
  term
}

##' Control the sequence of optimizations and termination criteria.
##'
##' @title Control parameters for optimization
##' @param grtol maximum gradient tolerance
##' @param dparrtol relative tolerance of parameter change between
##'   optimization calls
##' @param dobjrtol relative tolerance of objective value change between
##'   optimization calls
##' @param maxopts maximum number of optimizations to run
##' @param bfgs_control a list of control parameters for the BFGS optimizations
##' @return a \code{spatq_optcontrol} object
##' @author John K Best
##' @export
spatq_optcontrol <- function(grtol = 1e-8,
                             dparrtol = -Inf,
                             dobjrtol = -Inf,
                             maxopts = 1,
                             bfgs_control = list()) {
  structure(list(grtol = grtol,
                 dparrtol = dparrtol,
                 dobjrtol = dobjrtol,
                 maxopts = maxopts,
                 bfgs_control = bfgs_control),
            class = "spatq_optcontrol")
}

##' @describeIn fit_spatq Get object report
##' @export
report_spatq <- function(obj) {
  obj$report()
}

##' @describeIn fit_spatq Get object sdreport
##' @export
sdreport_spatq <- function(obj,
                           ## Don't bias correct if no random effects
                           bias.correct = !is.null(obj$env$random),
                           getJointPrecision = FALSE,
                           ...) {
  TMB::sdreport(obj,
                bias.correct = bias.correct,
                getJointPrecision = getJointPrecision,
                ...)
}

##' @describeIn fit_spatq Get finite difference Hessian
##' @importFrom stats optimHess
##' @export
hessian_spatq <- function(obj, fit) {
  optimHess(fit$par, obj$fn, obj$gr)
}
