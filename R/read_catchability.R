##' Read prep file and calculate realized catchability. Uses rocky habitat
##' covariate for \code{habq} and \code{bycatch} studies and
##' \code{log_catchability_devs} for other studies.
##'
##' @title Read catchability
##' @param study Study name
##' @param repl Replicate number
##' @param opmod Operating model number
##' @param root_dir Root directory
##' @param base_q Base catchability
##' @return List with \code{surv} and \code{comm} matrices
##' @author John K Best
##' @export
read_catchability <- function(study,
                              repl,
                              opmod,
                              root_dir = ".",
                              base_q = 0.01) {
  flnm <- sim_file_paths(study, repl, opmod, root_dir)$prep
  h5 <- hdf5r::h5file(flnm, mode = "r")
  if (study == "habq") {
    rocky <- h5[["rocky_habitat"]]$read()
    comm_q <- base_q * ifelse(rocky, 0.9, 1.0)
    surv_q <- base_q * ifelse(rocky, 0.1, 1.0)
  } else if (study == "bycatch") {
    rocky <- h5[["rocky_habitat"]]$read()
    comm_q <- base_q * ifelse(rocky, 0.5, 1.0)
    surv_q <- matrix(base_q, nrow = nrow(comm_q), ncol = ncol(comm_q))
  } else {
    logq <- h5[["log_catchability_devs"]][, , repl]
    ## Get logq scaling
    if (study == "qdevscaling") {
      s <- 10 ^ seq(-3, -0.5, 0.5)[opmod]
    } else {
      s <- 0.05
    }
    comm_q <- exp(s * logq - s ^ 2 / 2)
    surv_q <- matrix(base_q, nrow = nrow(comm_q), ncol = ncol(comm_q))
  }
  list(surv = surv_q,
       comm = comm_q)
}
