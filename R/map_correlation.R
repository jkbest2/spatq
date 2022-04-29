##' Calculates the correlation of true and estimated populations on a regular,
##' matching grid. It is expected that these are arrays with the first two
##' dimensions representing space and the third represents time. The
##' \code{byyear} argument allows for a correlation over the entire result, or
##' separate correlations for each year.
##'
##' @title Calculate the correlation between true and estimated populations
##' @param true Array of true biomass
##' @param est Array of estimated biomass, with matching first and second
##'   dimensions
##' @param byyear Calculate correlations by year?
##' @param method Which correlation coefficient to use; passed to
##'   \code{\link[stats]{cor}}
##' @return Correlation, or vector of correlations by year
##' @author John K Best
##' @export
map_correlation <- function(true, est, byyear = FALSE, method = "pearson") {
  all(dim(true)[1:2] == dim(est)[1:2]) || error("Domains must be the same size")
  nyears <- min(dim(true)[3], dim(est)[3])
  if (byyear) {
    mapcor <- vapply(1:nyears,
                     function(yr)
                       cor(c(true[, , yr]), c(est[, , yr]), method = method),
                     1.0)
  } else {
    mapcor <- cor(true[, , 1:nyears], est[, , 1:nyears], method = method)
  }
  return(mapcor)
}
