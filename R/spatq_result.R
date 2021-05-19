spatq_result <- function(studyspec, fit, lpb, rep, sdr) {
  structure(list(spec = studyspec,
                 fit = fit,
                 lpb = lpb,
                 rep = rep,
                 sdr = sdr),
            class = "spatq_result")
}

##' @title Maximum gradient component
##' @param x \code{spatq_fit} or \code{spatq_result}
##' @param abs Use absolute value of gradient?
##' @return Maximum gradient component of the fixed effect parameters
##' @author John K Best
mgc <- function(x, abs = TRUE) UseMethod("mgc")
##' @export
mgc.spatq_fit <- function(x, abs = TRUE) {
  grad <- x$grad
  if (abs) grad <- abs(grad)
  idx <- which.max(grad)
  gr <- grad[idx]
  names(gr) <- names(x$par)[idx]
  gr
}
##' @export
mgc.spatq_result <- function(x, abs = TRUE) {
  mgc(x$fit, abs)
}

##' @title Is the Hessian positive definite?
##' @param x \code{spatq_result} or \code{spatq_simstudyspec}
##' @return Logical
##' @author John K Best
##' @export
posdefhess <- function(x) UseMethod("posdefhess")
##' @export
posdefhess.spatq_simstudyspec <- function(x) {
  res <- read_rdata(x)
  posdefhess(res)
}
##' @export
posdefhess.spatq_result <- function(x) {
  x$sdr$pdHess
}
