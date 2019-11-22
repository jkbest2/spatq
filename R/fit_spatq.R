##' Use \code{optim} to estimate parameters of a \code{spatq} model.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param ... Additional options to pass to \code{optim}, e.g. \code{control}
##' @param method Optimization method to uses (see \code{optim} for available
##'   options)
##' @return An optimization object, report list, or sdreport list
##' @author John Best
##' @export
fit_spatq <- function(obj, ..., method = "BFGS") {
  stats::optim(obj$par, obj$fn, obj$gr, ...,
               method = method)
}

##' @describeIn fit_spatq Get object report
report_spatq <- function(obj) {
  obj$report()
}

##' @describeIn fit_spatq Get object sdreport
sdreport_spatq <- function(obj) {
  sdreport(obj)
}
