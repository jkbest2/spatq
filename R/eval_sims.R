##' Estimate the relationship between estimated indices and true biomass for
##' simulated data sets. Uses linear regression to estimate their relationship,
##' where one indicates unbiasedness. Requires a data frame with columns
##' 'est_index', 'true_biomass', and 'repl'.
##'
##' @title Estimate bias of estimation models
##' @param index_df Data frame with estimated indices, true biomass, and replicate
##'   number
##' @param unbiased Use debiased index estimates?
##' @return A list with 'alpha', the offsets for each replicate; 'delta', the
##'   bias coefficient; 'epsilon', the residuals; and 'sigma', the standard
##'   deviation of the regression.
##' @author John Best
##' @importFrom stats formula
##' @export
bias_metric <- function(index_df, unbiased = TRUE) {
  if (unbiased) {
    resp <- "raw_unb"
    form <- formula(log(raw_unb) ~ 0 + factor(repl) + log(raw_true))
  } else {
    resp <- "raw_est"
    form <- formula(log(raw_est) ~ 0 + factor(repl) + log(raw_true))
  }
  if (!all(c(resp, "raw_true", "repl") %in% names(index_df)))
    stop("Data frame must have columns \'raw_unb\', \'raw_true\', and \'repl\'")
  mod <- stats::lm(form, data = index_df)
  mod_coef <- stats::coef(mod)
  list(alpha = mod_coef[!grepl("raw_true", names(mod_coef))],
       delta = mod_coef[grepl("raw_true", names(mod_coef))],
       epsilon = stats::resid(mod),
       sigma = stats::sigma(mod))

}

##' Rescale indices of abundance so that they have product one (mean-log zero).
##' This facilitates comparison among models that only estimate relative
##' abundance and with simulated true abundances
##'
##' @title Rescale indices
##' @param b Unscaled indices/biomass
##' @param sd Unscaled index standard deviation
##' @return A list with elements \code{index} and \code{sd}, where indices are
##'   scaled to have a product 1 (mean log of zero) and standard deviations are
##'   scaled by the same factor.
##' @author John Best
##' @export
rescale_index <- function(b, sd = NULL) {
  scale <- exp(mean(log(b)))
  index <- list(index = b / scale)
  if (!is.null(sd)) {
    sc_sd <- sd / scale
    index$sd <- sc_sd
  }
  attr(index, "scale") <- scale
  class(index) <- "spatq_index"
  index
}

##' Calculate the root mean square error of the estimated indices of abundance
##' compared to the true, simulated values. Expects that these will be on a
##' common scale, e.g. by using \code{rescale_indices}.
##'
##' @title Root mean square error
##' @param est_index Estimated rescaled index
##' @param true_index True rescaled index
##' @return Root mean square error of the estimates
##' @author John Best
##' @export
rmse_metric <- function(est_index, true_index) {
  sqrt(mean((est_index - true_index)^2))
}

##' In order to assess calibration of the estimation models, find the percentile
##' of the true index under the estimated distribution, assumed to be normal
##' with estimated mean and standard deviation.
##'
##' @title Percentile of true value in estimated distribution
##' @param true_index True rescaled index
##' @param est_index Estimated index of abundance
##' @param est_sd Estimated standard deviation of index
##' @return The percentile of the true value
##' @author John Best
##' @export
true_percentile <- function(true_index, est_index, est_sd) {
  stats::pnorm(true_index, est_index, est_sd)
}
