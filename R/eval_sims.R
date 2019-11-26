##' Estimate the relationship between estimated indices and true biomass for
##' simulated data sets. Uses linear regression to estimate their relationship,
##' where one indicates unbiasedness. Requires a data frame with columns
##' 'est_index', 'true_biomass', and 'repl'.
##'
##' .. content for \details{} ..
##' @title Estimate bias of estimation models
##' @param fit_df Data frame with estimated indices, true biomass, and replicate
##'   number
##' @return A list with 'alpha', the offsets for each replicate; 'delta', the
##'   bias coefficient; 'epsilon', the residuals; and 'sigma', the standard
##'   deviation of the regression.
##' @author John Best
##' @export
bias_metric <- function(fit_df) {
  if (any(c("est_index", "true_biomass", "repl") %ni% names(fit_df)))
    stop("Data frame must have columns \'est_index\', \'true_biomass\', and \'repl\'")
  mod <- lm(log(est_index) ~ 0 + factor(repl) + log(true_biomass),
            data = fit_df)
  mod_coef <- coef(mod)
  list(alpha = mod_coef[!grepl("true_biomass", names(mod_coef))],
       delta = mod_coef[grepl("true_biomass", names(mod_coef))],
       epsilon = resid(mod),
       sigma = sigma(mod))

}

##' Rescale indices of abundance so that they have product one (mean-log zero).
##' This facilitates comparison among models that only estimate relative
##' abundance and with simulated true abundances
##'
##' @title Rescale indices
##' @param b Unscaled indices/biomass
##' @return A vector of indices scaled to have a product 1 (mean log of zero)
##' @author John Best
##' @export
rescale_index <- function(b) {
  scale <- exp(mean(log(b)))
  i <- b / i
  attr(i, "scale") <- scale
  i
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
  pnorm(true_index, est_index, est_sd)
}
