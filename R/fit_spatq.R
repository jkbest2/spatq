##' Use \code{optim} to estimate parameters of a \code{spatq} model.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param ... Additional options to pass to \code{optim}, e.g. \code{control}
##' @return An optimization object, report list, or sdreport list
##' @author John Best
##' @export
fit_spatq <- function(obj, ...) {
  stats::optim(obj$par, obj$fn, obj$gr, ...)
}

##' @describeIn fit_spatq Get object report
report_spatq <- function(obj) {
  obj$report()
}

##' @describeIn fit_spatq Get object sdreport
sdreport_spatq <- function(obj) {
  sdreport(obj)
}
