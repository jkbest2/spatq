##' Use \code{optimx::Rcgmin} to find quickly improve, then
##' \code{optimx::Rvmmin} to improve once closer to optimum. After each 50
##' iterations of \code{optimx::Rvmmin} check the maximum gradient component,
##' stopping if it is less than \code{grtol}.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param fit Previous fit to use as starting values
##' @param ... Additional options to pass to \code{optim}, e.g. \code{control}
##' @param grtol Gradient tolerance
##' @param maxopts Maximum number of variable metric optimizations to run. 50
##'   iterations each and checking maximum gradient component.
##' @return An optimization object, report list, or sdreport list
##' @author John Best
##' @export
fit_spatq <- function(obj, fit = NULL, ..., grtol = 1e-3, maxopts = 10) {
  ## Use default starting parameters if previous fit is not provided
  if (is.null(fit)) {
    fit <- optimx::Rcgmin(obj$par, obj$fn, obj$gr, control = list(maxit = 50))
  }
  for (i in 1:maxopts) {
    fit <- optimx::Rvmmin(fit$par, obj$fn, obj$gr, ...)
    attr(fit, "mgc") <- max(obj$gr(fit$par))
    if (attr(fit, "mgc") < grtol) break
  }
  fit
}

##' @describeIn fit_spatq Get object report
##' @export
report_spatq <- function(obj) {
  obj$report()
}

##' @describeIn fit_spatq Get object sdreport
##' @export
sdreport_spatq <- function(obj) {
  TMB::sdreport(obj)
}
