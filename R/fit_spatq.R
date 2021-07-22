##' Use nonlinear optimization routines from \code{\link[stats]{nlminb}} or
##' \code{\link[stats]{optim}} to fit. Can include constraints where appropriate
##' for the chosen method. If \code{method = "nlminb"} is used, the result is
##' restructure to match \code{\link[stats]{optim}} output using
##' \code{\link{fix_nlminb_fit}}.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A \code{\link{spatq_obj}}
##' @param fit Previous fit, use results as starting values
##' @param method Optimization method to use, as a string. Valid options are
##'   \code{"nlminb"} or one of the multivariate \code{\link[stats]{optim}}
##'   methods.
##' @param bounds List with \code{"upper"} and \code{"lower"} parameter bounds;
##'   will error for unconstrained optimization methods
##' @param control list of control arguments to pass to
##'   \code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}
##' @return A spatq_fit object containing the optimization output and
##'   optimization diagnostics
##' @importFrom stats optim nlminb
##' @author John K Best
##' @export
spatq_fit <- function(obj,
                      fit = NULL,
                      method = "nlminb",
                      bounds = list(lower = -Inf, upper = Inf),
                      control = list()) {
  if (is.null(fit)) {
    fit <- init_spatq_fit(obj)
  }
  t0 <- Sys.time()
  if (method == "nlminb") {
    newfit <- nlminb(fit$par, obj$fn, obj$gr,
                  control = control)
    ## Change resulting fit to match `optim` output
    newfit <- fix_nlminb_fit(newfit)
  } else if (method  %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) {
    newfit <- optim(fit$par, obj$fn, obj$gr,
                    method = method,
                    lower = bounds$lower, upper = bounds$upper,
                    control = control)
  } else {
    stop("Method ", method, " not available, use nlminb or optim method")
  }
  newfit$opt_time <- Sys.time() - t0
  fit <- attach_optdiags(newfit, fit, obj)
  return(fit)
}

##' Reorders output, renames "objective" to "value", and drops "iterations".
##'
##' @title Make \code{nlminb} output match \code{optim} output
##' @param fit \code{\link[stats]{nlminb}} output
##' @return Optimization result that matches \code{\link[stats]{optim}}
##' @author John K Best
##' @export
fix_nlminb_fit <- function(fit) {
  list(par = fit$par,
       value = fit$objective,
       counts = fit$evaluations,
       convergence = fit$convergence,
       message = fit$message)
}

##' Uses \code{\link[stats]{optim}} to find an optimum. If a previous fit is not
##' provided, the first fit uses the conjugate gradient method to improve
##' quickly. After that, the BFGS method is used to refine the optimization.
##' Termination tolerances and control parameters for \code{optim} can be
##' provided through \code{\link{spatq_optcontrol}}.
##'
##' @title Fit a spatq model by maximum likelihood
##' @param obj A spatq \code{ADFun}, as returned by \code{prepare_adfun} or
##'   \code{make_sim_adfun}
##' @param fit Previous fit to use as starting values
##' @param optcontrol a
##' @param bias.correct Use bias correction for \code{\link[TMB]{sdreport}}?
##' @param bias.correct.control Control list for \code{bias.correct}
##' @param getJointPrecision Return the joint fixed and random effect precision
##'   matrix from \code{\link[TMB]{sdreport}}?
##' @param ... additional arguments to pass to \code{\link[TMB]{sdreport}}
##' @return An optimization object, report list, or sdreport list
##' @author John Best
##' @importFrom stats optim
##' @export
fit_spatq <- function(obj, fit = NULL, optcontrol = spatq_optcontrol()) {
  UseMethod("fit_spatq")
}
##' @export
fit_spatq.spatq_obj <- function(obj, fit = NULL, optcontrol = spatq_optcontrol()) {
  if (is.null(fit)) {
    fit <- init_spatq_fit(obj)
  }
  while (!terminate_opt(fit, optcontrol)) {
    fit1 <- optim(fit$par, obj$fn, obj$gr, method = "BFGS",
                  control = optcontrol$bfgs_control)
    fit <- attach_optdiags(fit1, fit, obj)
  }
  fit
}
##' @export
fit_spatq.spatq_designobj <- function(obj, fit = NULL, optcontrol = spatq_optcontrol()) {
  index_df <- obj %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(index = mean(catch_biomass),
                     n = dplyr::n(),
                     sd = sd(catch_biomass) / sqrt(n))
  index <- index_df$index
  names(index) <- purrr::rep_along(index, "Index")
  sd <- index_df$sd

  structure(list(value = index,
                 sd = sd,
                 unbiased = list(value = index,
                                 sd = sd)),
            class = "spatq_designfit")
}

##' Attach optimization diagnostics comparing a previous optimization result to
##' a new optimization result.
##'
##' @title Attach optimization diagnostics
##' @param newfit new optimization result
##' @param oldfit previous optimization result
##' @param obj TMB objective function
##' @return A \code{spatq_fit} object
##' @author John K Best
##' @export
attach_optdiags <- function(newfit, oldfit = NULL, obj) {
  if (is.null(oldfit)) {
    oldfit <- init_spatq_fit(obj)
  }
  newfit$grad <- obj$gr(newfit$par)
  newfit$dgrad <- oldfit$grad - newfit$grad
  newfit$dparrel <- (oldfit$par - newfit$par) / oldfit$par
  newfit$dobjrel <- (oldfit$value - newfit$value) / oldfit$value
  newfit$nfit <- oldfit$nfit + 1
  newfit$totcounts <- oldfit$totcounts + newfit$counts
  newfit$opt_time <- oldfit$opt_time + newfit$opt_time
  class(newfit) <- "spatq_fit"
  newfit
}

##' @method print spatq_fit
##' @author John K Best
##' @export
print.spatq_fit <- function(x, ...) {
  fit_df <- data.frame(parname = names(x$par),
                       parval = x$par,
                       grad = c(x$grad))
  print(fit_df)
  cat("\nObjective: ", x$value, "\n",
      "Convergence: ", x$convergence, "\n",
      "Message: ", x$message, "\n",
      "MGC: ", max(abs(x$grad)), "\n",
      "Time optimizing: ", x$opt_time, units(x$opt_time), "\n",
      sep = "")

  return(invisible(x))
}

##' @describeIn attach_optdiags Initialize a spatq_fit object
##' @export
init_spatq_fit <- function(obj) {
  ## Fill in with values that make sense or else values that won't cause
  ## optimization to terminate prematurely
  structure(
    list(value = obj$fn(obj$par),
         par = obj$par,
         grad = obj$gr(obj$par),
         dparrel = rep_len(Inf, length(obj$par)),
         dobjrel = Inf,
         nfit = 0,
         totcounts = c("function" = 0, "gradient" = 0),
         opt_time = as.difftime("0", "%S")),
    class = "spatq_fit")
}

##' Sequentially tests whether maximum number of optimizations has been reached,
##' if the maximum gradient component is less than \code{optcontrol$grtol},
##' relative parameter value changes between optimizations is less than
##' \code{optcontrol$dparrtol}, or relative change in the objective function
##' value is less than \code{optcontrol$dobjrtol}.
##'
##' @title Test whether optimization should be stopped
##' @param fit \code{spatq_fit} object from \code{\link{attach_optdiags}}
##' @param optcontrol \code{\link{spatq_optcontrol}} object
##' @return boolean indicating termination (TRUE) or not (FALSE)
##' @author John K Best
##' @export
terminate_opt <- function(fit = NULL, optcontrol) {
  term <- FALSE
  if (is.null(fit)) {
    term <- FALSE
  } else if (fit$nfit >= optcontrol$maxopts) {
    term <- TRUE
  ## } else if (max(abs(fit$grad)) < optcontrol$grtol) {
  ##   term <- TRUE
  ## } else if (max(abs(fit$dparrel)) < optcontrol$dparrtol) {
  ##   term <- TRUE
  ## } else if (fit$dobjrel < optcontrol$dobjrtol) {
  ##   term <- TRUE
  }
  term
}

##' Control the sequence of optimizations and termination criteria.
##'
##' @title Control parameters for optimization
##' @param grtol maximum gradient tolerance
##' @param dparrtol relative tolerance of parameter change between
##'   optimization calls
##' @param dobjrtol relative tolerance of objective value change between
##'   optimization calls
##' @param maxopts maximum number of optimizations to run
##' @param bfgs_control a list of control parameters for the BFGS optimizations
##' @return a \code{spatq_optcontrol} object
##' @author John K Best
##' @export
spatq_optcontrol <- function(grtol = 1e-8,
                             dparrtol = -Inf,
                             dobjrtol = -Inf,
                             maxopts = 1,
                             bfgs_control = list()) {
  structure(list(grtol = grtol,
                 dparrtol = dparrtol,
                 dobjrtol = dobjrtol,
                 maxopts = maxopts,
                 bfgs_control = bfgs_control),
            class = "spatq_optcontrol")
}

##' @describeIn fit_spatq Get object report
##' @export
report_spatq <- function(obj) {
  UseMethod("report_spatq")
}
##' @export
report_spatq.spatq_obj <- function(obj) {
  obj$report()
}
##' @export
report_spatq.spatq_designobj <- function(obj) {
  list()
}

##' @describeIn fit_spatq Get object sdreport
##' @export
sdreport_spatq <- function(obj,
                           ## Don't bias correct if no random effects
                           bias.correct = !is.null(obj$env$random),
                           bias.correct.control = list(sd = TRUE),
                           getJointPrecision = FALSE,
                           ...) {
  UseMethod("sdreport_spatq")
}
##' @export
sdreport_spatq.spatq_obj <- function(obj,
                                     ## Don't bias correct if no random effects
                                     bias.correct = !is.null(obj$env$random),
                                     bias.correct.control = list(sd = TRUE),
                                     getJointPrecision = FALSE,
                                     ...) {
  TMB::sdreport(obj,
                bias.correct = bias.correct,
                bias.correct.control = bias.correct.control,
                getJointPrecision = getJointPrecision,
                ...)
}
##' @export
sdreport_spatq.spatq_designobj <- function(obj) {
  index <- fit_spatq(obj)
  ## class(index) <- c("sdreport_spatq")
  index
}

##' @describeIn fit_spatq Get finite difference Hessian
##' @importFrom stats optimHess
##' @export
hessian_spatq <- function(obj, fit) {
  optimHess(fit$par, obj$fn, obj$gr)
}
