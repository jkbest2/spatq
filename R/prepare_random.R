##' Prepare a character vector indicating which parameters should be
##' marginalized out using the Laplace approximation. Does not include
##' parameters that are \code{map}'d.
##'
##' Parameters in this model that are typically marginalized are:
##' \itemize{
##'   \item \code{gamma_n}, \code{gamma_w}
##'   \item \code{omega_n}, \code{omega_w}
##'   \item \code{epsilon_n}, \code{epsilon_w}
##'   \item \code{eta_n}, \code{eta_w}
##'   \item \code{phi_n}, \code{phi_w}
##'   \item \code{psi_n}, \code{psi_w}
##' }
##'
##' @title Prepare \code{random}
##' @param map A \code{map} list, as from \code{prepare_map}
##' @return A character vector indicating which parameters should be integrated
##'   out using the Laplace approximation.
##' @author John Best
##' @export
prepare_random <- function(map) {
  re_pars <- c("gamma_n", "gamma_w",
               "omega_n", "omega_w",
               "epsilon_n", "epsilon_w",
               "eta_n", "eta_w",
               "phi_n", "phi_w",
               "psi_n", "psi_w")
  ## Check that all `map` random effects parameter entries are all NAs; check
  ## that none are included but not actually map'd
  if (!all(vapply(re_pars, function(p) all(is.na(map[[p]])),
                  FUN.VALUE = TRUE)))
    stop("Map'd random effects parameters must be all NAs")
  ## Return vector of random effects names with map'd parameter names removed
  setdiff(re_pars, names(map))
}
