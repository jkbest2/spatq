% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_spatq.R
\name{spatq_optcontrol}
\alias{spatq_optcontrol}
\title{Control parameters for optimization}
\usage{
spatq_optcontrol(
  grtol = 1e-08,
  dparrtol = -Inf,
  dobjrtol = -Inf,
  maxopts = 1,
  bfgs_control = list()
)
}
\arguments{
\item{grtol}{maximum gradient tolerance}

\item{dparrtol}{relative tolerance of parameter change between
optimization calls}

\item{dobjrtol}{relative tolerance of objective value change between
optimization calls}

\item{maxopts}{maximum number of optimizations to run}

\item{bfgs_control}{a list of control parameters for the BFGS optimizations}
}
\value{
a \code{spatq_optcontrol} object
}
\description{
Control the sequence of optimizations and termination criteria.
}
\author{
John K Best
}
