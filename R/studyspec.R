##' Specification of a single simulation study fit.
##'
##' @title Simulation study specification
##' @param s list of spec components (arguments listed below) or study name
##' @param ... Rest of the spec components if s is study name
##' @param study Simulation study name
##' @param repl Replicate number
##' @param opmod Operating model number
##' @param estmod Estimation model name
##' @param sub_df Data subsetting data frame
##' @param estd Parameter estimation specification
##' @param root_dir Root directory
##' @return The same list, but with class \code{spatq_simstudyspec}. Only checks
##'   that all elements are non-\code{NULL}.
##' @author John K Best
##' @export
spatq_simstudyspec <- function(s, ...) UseMethod("spatq_simstudyspec")
##' @export
spatq_simstudyspec.list <- function(s, ...) {
  new_spatq_simstudyspec(s$study,
                         s$repl,
                         s$opmod,
                         s$estmod,
                         s$sub_df,
                         s$estd,
                         s$root_dir)
}
##' @export
spatq_simstudyspec.data.frame <- function(s, ...) {
  spatq_simstudyspec(as.list(s))
}

##' @export
spatq_simstudyspec.character <- function(s, ...) {
  new_spatq_simstudyspec(s, ...)
}

##' @describeIn spatq_simstudyspec New simstudyspec
##' @export
new_spatq_simstudyspec <- function(study, repl, opmod, estmod, sub_df = NULL, estd = NULL, root_dir = ".") {
  ## Set default root as current working directory
  if (is.null(root_dir)) root_dir <- "."
  spec <- structure(list(study = study,
                         repl = repl,
                         opmod = opmod,
                         estmod = estmod,
                         sub_df = sub_df,
                         estd = estd,
                         root_dir = root_dir),
                    class = "spatq_simstudyspec")
  if (any(vapply(spec[c(1:4, 7)], is.null, TRUE))) {
    stop("simstudyspec is missing pieces")
  }
  spec
}
