% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studyspec.R
\name{spatq_simstudyspec}
\alias{spatq_simstudyspec}
\alias{new_spatq_simstudyspec}
\title{Simulation study specification}
\usage{
spatq_simstudyspec(s, ...)

new_spatq_simstudyspec(
  study,
  repl,
  opmod,
  estmod,
  sub_df = NULL,
  estd = NULL,
  root_dir = "."
)
}
\arguments{
\item{s}{list of spec components (arguments listed below) or study name}

\item{...}{Rest of the spec components if s is study name}

\item{study}{Simulation study name}

\item{repl}{Replicate number}

\item{opmod}{Operating model number}

\item{estmod}{Estimation model name}

\item{sub_df}{Data subsetting data frame}

\item{estd}{Parameter estimation specification}

\item{root_dir}{Root directory}
}
\value{
The same list, but with class \code{spatq_simstudyspec}. Only checks
that all elements are non-\code{NULL}.
}
\description{
Specification of a single simulation study fit.
}
\section{Functions}{
\itemize{
\item \code{new_spatq_simstudyspec}: New simstudyspec
}}

\author{
John K Best
}
