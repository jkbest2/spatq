#' Spatially varying catchability for fishery-dependent observations
#'
#' A package for fitting index standardization models that combine
#' fishery-independent and fishery-dependent observations, where the
#' fishery-dependent observations are allowed a spatially varying
#' catchability.
#'
#' @docType package
#' @name spatq
#' @useDynLib spatq
NULL

## These variable names are used in NSE functions and throw R CMD check warnings
## if they're not declared here. They could probably be converted to standard
## evaluation using `rlang`, but that's for later.
globalVariables(c(
  "catch_biomass",
  "cloglog",
  "coord1",
  "coord2",
  "data",
  "estmod",
  "s1",
  "loc_idx",
  "n",
  "opmod",
  "pop",
  "raw_est",
  "raw_sd",
  "raw_true",
  "raw_unb",
  "raw_unb_sd",
  "repl",
  "sample_n",
  "time",
  "vessel_idx",
  "year"))
