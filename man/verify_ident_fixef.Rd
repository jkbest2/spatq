% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verify_spatq.R
\name{verify_ident_fixef}
\alias{verify_ident_fixef}
\title{Verify that fixed effects are identifiable}
\usage{
verify_ident_fixef(data)
}
\arguments{
\item{data}{Data list, as produced by \code{prepare_data}}
}
\value{
TRUE if it passes, error otherwise
}
\description{
Check that the combined fixed effect design matrices for abundance and
catchability are full rank.
}
\author{
John Best
}
