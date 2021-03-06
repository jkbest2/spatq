% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_fits.R
\name{save_fit}
\alias{save_fit}
\title{Save fit results}
\usage{
save_fit(studyspec, fit, lpb, rep, sdr, root_dir = ".")
}
\arguments{
\item{studyspec}{A \code{spatq_studyspec} specification}

\item{fit}{A \code{\link{spatq_fit}} fit result}

\item{lpb}{The contents of \code{obj$env$last.par.best}, gathered into a
named list}

\item{rep}{A parameter report from \code{\link{report_spatq}}}

\item{sdr}{An SD Report from \code{\link{sdreport_spatq}}}

\item{root_dir}{Root directory}
}
\value{
The name of the saved file
}
\description{
Save the results of a fit.
}
\author{
John K Best
}
