##' Specification of a single simulation study fit.
##'
##' @title Simulation study specification
##' @param speclist List with \code{study}, \code{repl}, \code{opmod},
##'   \code{estmod}, \code{sub_df}, \code{estd}
##' @param study Simulation study name
##' @param repl Replicate number
##' @param opmod Operating model number
##' @param estmod Estimation model name
##' @param sub_df Data subsetting data frame
##' @param estd Parameter estimation specification
##' @return The same list, but with class \code{spatq_simstudyspec}. Only checks
##'   that all elements are non-\code{NULL}.
##' @author John K Best
##' @export
spatq_simstudyspec <- function(speclist, feather = TRUE) {
  new_spatq_simstudyspec(speclist$study,
                         speclist$repl,
                         speclist$opmod,
                         speclist$estmod,
                         speclist$sub_df,
                         speclist$estd,
                         speclist$root_dir)
}

##' @describeIn spatq_simstudyspec New simstudyspec
##' @export
new_spatq_simstudyspec <- function(study, repl, opmod, estmod, sub_df, estd, root_dir) {
  spec <- structure(list(study = study,
                         repl = repl,
                         opmod = opmod,
                         estmod = estmod,
                         sub_df = sub_df,
                         estd = estd,
                         root_dir = root_dir),
                    class = "spatq_simstudyspec")
  if (any(vapply(spec, is.null, TRUE))) {
    stop("simstudyspec is missing pieces")
  }
  spec
}
