##' Return a data frame compatible with \code{\link{subsample_catch}} specifying
##' the number of catch observations to use for each estimation model.
##'
##' Valid estimation model names that subset observations to use only survey data are:
##'
##' \itemize{
##'   \item \code{"design"}
##'   \item \code{"model"}
##'   \item \code{"survey"}
##'   \item \code{"survey_spt"}
##' }
##'
##' Estimation models that include all survey observations and \code{n_comm} commercial observations
##' are:
##' \itemize{
##'   \item \code{"design_all"}
##'   \item \code{"model_all"}
##'   \item \code{"spatial_ab"}
##'   \item \code{"spatial_q"}
##'   \item \code{"sptemp_ab"}
##'   \item \code{"sptemp_q"}
##' }
##' @title Specify data subsample for estimation model
##' @param estmod Estimation model name
##' @param n_comm Number of commercial observations to use
##' @return A data frame to be used by \code{\link{subsample_catch}}
##' @author John K Best
##' @export
em_subsample <- function(estmod, n_comm = 4000) UseMethod("em_subsample")
##' @export
em_subsample.character <- function(estmod, n_comm = 4000) {
  switch(estmod,
         design = ,
         model = ,
         survey = ,
         survey_spt = data.frame(vessel_idx = 2, n = 0),
         design_all = ,
         model_all = ,
         spatial_ab = ,
         spatial_q = ,
         sptemp_ab = ,
         sptemp_q = data.frame(vessel_idx = 2, n = n_comm),
         stop("Invalid estimation model"))
}
##' @export
em_subsample.data.frame <- function(estmod, n_comm = 4000) estmod
##' @export
em_subsample.NULL <- function(estmod, n_comm = 4000) NULL
##' @export
em_subsample.default <- function(estmod, n_comm = 4000) {
  stop("Must pass either data frame or name of valid estimation model")
}

##' Estimated parameter list for a set of common estimation models. Otherwise
##' will pass through a specification.
##'
##' Valid estimation models are:
##'
##' \describe{
##'   \item{\code{design}}{Only yearly intercepts}
##'   \item{\code{model}}{Only yearly intercepts}
##'   \item{\code{survey}}{Yearly intercepts and spatial abundance}
##'   \item{\code{survey_spt}}{Yearly intercepts with spatial and spatiotemporal abundance}
##'   \item{\code{design_all}}{Yearly intercepts}
##'   \item{\code{model_all}}{Yearly and vessel intercepts}
##'   \item{\code{spatial_ab}}{Yearly and vessel intercepts, spatial abundance}
##'   \item{\code{spatial_q}}{Yearly and vessel intercepts, spatial abundance and catchability}
##'   \item{\code{sptemp_ab}}{Yearly and vessel intercepts with spatial and spatiotemporal abundance}
##'   \item{\code{sptemp_q}}{Yearly and vessel intercepts, spatiotemporal abundance and spatial catchability}
##' }
##'
##' Default is to use a Tweedie observation likelihood and to map all
##' \code{kappa} decorrelation rate parameters to the same value.
##' @title Estimated parameters for an estimation model
##' @param estmod Name of an estimation model or a \code{spatq_spec}
##'   object
##' @param obs_lik Use delta Poisson link log-normal (\code{0}) or Tweedie
##'   (\code{1})?
##' @return A \code{spatq_spec} object
##' @author John K Best
##' @export
em_estd <- function(estmod, obs_lik = 1L) UseMethod("em_estd")
##' @export
em_estd.spatq_spec <- function(estmod, obs_lik = 1L) estmod
##' @export
em_estd.character <- function(estmod, obs_lik = 1L) {
  switch(estmod,
         design = ,
         model = specify_estimated(beta = TRUE,
                                   lambda = FALSE,
                                   obs_lik = obs_lik),
         survey = specify_estimated(beta = TRUE,
                                    omega = TRUE,
                                    lambda = FALSE,
                                    obs_lik = obs_lik),
         survey_spt = specify_estimated(beta = TRUE,
                                        omega = TRUE,
                                        epsilon = TRUE,
                                        lambda = FALSE,
                                        kappa_map = c(1, NA,
                                                      1, NA,
                                                      NA, NA,
                                                      NA, NA),
                                        obs_lik = obs_lik),
         design_all = specify_estimated(beta = TRUE,
                                        lambda = FALSE,
                                        obs_lik = obs_lik),
         model_all = specify_estimated(beta = TRUE,
                                       lambda = TRUE,
                                       obs_lik = obs_lik),
         spatial_ab = specify_estimated(beta = TRUE,
                                        omega = TRUE,
                                        lambda = TRUE,
                                        obs_lik = obs_lik),
         spatial_q = specify_estimated(beta = TRUE,
                                       omega = TRUE,
                                       lambda = TRUE,
                                       phi = TRUE,
                                       kappa_map = c(1, NA,
                                                     NA, NA,
                                                     1, NA,
                                                     NA, NA),
                                       obs_lik = obs_lik),
         sptemp_ab = specify_estimated(beta = TRUE,
                                       omega = TRUE,
                                       epsilon = TRUE,
                                       lambda = TRUE,
                                       kappa_map = c(1, NA,
                                                     1, NA,
                                                     NA, NA,
                                                     NA, NA),
                                       obs_lik = obs_lik),
         sptemp_q = specify_estimated(beta = TRUE,
                                      omega = TRUE,
                                      epsilon = TRUE,
                                      lambda = TRUE,
                                      phi = TRUE,
                                      kappa_map = c(1, NA,
                                                    1, NA,
                                                    1, NA,
                                                    NA, NA),
                                      obs_lik = obs_lik),
         stop("Invalid estimation model"))
}
##' @export
em_estd.default <- function(estmod, obs_lik = 1L) {
  stop("Must pass either spatq_spec or name of valid estimation model")
}
