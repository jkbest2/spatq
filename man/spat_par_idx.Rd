% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_pars.R
\name{spat_par_idx}
\alias{spat_par_idx}
\title{Spatial parameter index}
\usage{
spat_par_idx(par_name)
}
\arguments{
\item{par_name}{Name of the spatial/spatiotemporal parameter (e.g.
\code{omega_n}, \code{phi_w} &c.)}
}
\value{
The index of the effect in the parameter vectors
}
\description{
Spatial/spatiotemporal parameters (kappa and tau) are stored in a length-8
vector. This function takes a name and returns the corresponding index.
}
\author{
John Best
}
