##' Prepare a \code{map} argument for \code{MakeADFun}.
##'
##' @title Prepare a \code{map}
##' @param pars Parameter list, as from \code{prepare_pars}
##' @param spec List of logicals indicating which parameters are to be
##'   estimated, as output by \code{\link{specify_estimated}}
##' @return A \code{map} list suitable for passing to \code{MakeADFun}
##' @author John Best
##' @export
prepare_map <- function(pars, spec) {
  ## Drop parameters that are *not* mapped
  spec2 <- spec[!vapply(spec, all, TRUE)]
  ## Invert to specify *mapped* parameter vectors of correct length
  mapped <- lapply(
    names(spec2),
    function(nm) {
      rep(!spec2[[nm]], length.out = length(pars[[nm]]))
    }
  )
  names(mapped) <- names(spec2)
  map <- lapply(
    mapped,
    function(mpd) {
      ## Add individual levels for unmapped parameters
      v <- cumsum(!mpd)
      ## Specify NA for mapped parameters
      v[mpd] <- NA
      ## Convert to factor vector for TMB
      factor(v)
    }
  )

  ## If kappa map is explicitly provided, replace the default
  if (!is.null(attr(spec, "kappa_map"))) {
    kappa_map <- attr(spec, "kappa_map")
    if (any(xor(is.na(map$log_kappa), is.na(kappa_map)))) {
      stop("kappa_map must include NAs in correct locations")
    }
    map$log_kappa <- factor(kappa_map)
  }

  ## Anisotropy parameters default to identity matrix
  if (!spec$H_pars) {
    map$H_pars <- factor(c(NA, NA))
  }

  ## Drop `obs_lik`; currently always estimating these parameters
  map$obs_lik <- NULL

  return(map)
}
