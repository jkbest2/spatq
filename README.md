# spatq

<!-- badges: start -->
<!-- badges: end -->

This package provides a [Template Model Builder](https://github.com/kaskr/adcomp) program to fit index standardization models that include both fishery-dependent and -independent observations. The goal is to explore the utility of standardizing fishery-dependent observations by allowing for spatially varying catchability.

## Installation

This package has not been released. The most recent version can be installed with:

``` r
devtools::install_github("jkbest2/spatq")
```

## Example

A basic example of generating fake data and fitting a simple model is available in the `tests/testthat/test-spatq.R` file. Currently the `data` and `parameters` arguments to `MakeADFun` must be built by the user. This must be done carefully to avoid mismatched dimensions that can crash your R session.

