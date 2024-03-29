% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_specs.R
\name{em_estd}
\alias{em_estd}
\title{Estimated parameters for an estimation model}
\usage{
em_estd(estmod, obs_lik = 1L)
}
\arguments{
\item{estmod}{Name of an estimation model or a \code{spatq_spec}
object}

\item{obs_lik}{Use delta Poisson link log-normal (\code{0}) or Tweedie
(\code{1})?}
}
\value{
A \code{spatq_spec} object
}
\description{
Estimated parameter list for a set of common estimation models. Otherwise
will pass through a specification.
}
\details{
Valid estimation models are:

\describe{
\item{\code{design}}{Only yearly intercepts}
\item{\code{model}}{Only yearly intercepts}
\item{\code{survey}}{Yearly intercepts and spatial abundance}
\item{\code{survey_spt}}{Yearly intercepts with spatial and spatiotemporal abundance}
\item{\code{design_all}}{Yearly intercepts}
\item{\code{model_all}}{Yearly and vessel intercepts}
\item{\code{spatial_ab}}{Yearly and vessel intercepts, spatial abundance}
\item{\code{spatial_q}}{Yearly and vessel intercepts, spatial abundance and catchability}
\item{\code{sptemp_ab}}{Yearly and vessel intercepts with spatial and spatiotemporal abundance}
\item{\code{sptemp_q}}{Yearly and vessel intercepts, spatiotemporal abundance and spatial catchability}
}

Default is to use a Tweedie observation likelihood and to map all
\code{kappa} decorrelation rate parameters to the same value.
}
\author{
John K Best
}
