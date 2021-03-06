% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_data.R
\name{prepare_proc_switch}
\alias{prepare_proc_switch}
\title{Prepare process switch}
\usage{
prepare_proc_switch(random)
}
\arguments{
\item{random}{Character vector of random parameters; as from
\code{prepare_random}}
}
\value{
Logical vector of length 12 indicating which random processes are
\emph{not} map'd off
}
\description{
Generates a vector of length 6 indicating whether each pair of numbers
density and weight per group random processes should be included in the
likelihood calculation; i.e. they are map'd.
}
\author{
John Best
}
