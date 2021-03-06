% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{par_expr}
\alias{par_expr}
\title{Pretty expressions for plotting parameter names}
\usage{
par_expr(parname, hat = FALSE)
}
\arguments{
\item{parname}{Parameter name as a string}

\item{hat}{Hat for estimated parameter?}
}
\value{
An \code{expression} representing the parameter
}
\description{
Easy way to add parameters as Greek letters with subscripts to plots. Not
exported.
}
\author{
John K Best
}
