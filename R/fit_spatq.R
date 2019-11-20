##' Use \code{optim} to estimate parameters of a \code{spatq} model.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param ... Additional options to pass to \code{optim}, e.g. \code{control}
##' @return An optimiz
##' @author John Best
fit_spatq <- function(obj, ...) {
  optim(obj$par, obj$fn, obj$gr, ...)
}
