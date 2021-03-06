% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_spatq.R
\name{spatq_fit}
\alias{spatq_fit}
\title{Fit a spatq model by maximum likelihood}
\usage{
spatq_fit(
  obj,
  fit = NULL,
  method = "nlminb",
  bounds = list(lower = -Inf, upper = Inf),
  control = list()
)
}
\arguments{
\item{obj}{A \code{\link{spatq_obj}}}

\item{fit}{Previous fit, use results as starting values}

\item{method}{Optimization method to use, as a string. Valid options are
\code{"nlminb"} or one of the multivariate \code{\link[stats]{optim}}
methods.}

\item{bounds}{List with \code{"upper"} and \code{"lower"} parameter bounds;
will error for unconstrained optimization methods}

\item{control}{list of control arguments to pass to
\code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}}
}
\value{
A spatq_fit object containing the optimization output and
optimization diagnostics
}
\description{
Use nonlinear optimization routines from \code{\link[stats]{nlminb}} or
\code{\link[stats]{optim}} to fit. Can include constraints where appropriate
for the chosen method. If \code{method = "nlminb"} is used, the result is
restructure to match \code{\link[stats]{optim}} output using
\code{\link{fix_nlminb_fit}}.
}
\author{
John K Best
}
