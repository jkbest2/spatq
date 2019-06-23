# spatq

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jkbest2/spatq.svg?branch=master)](https://travis-ci.org/jkbest2/spatq)
[![Codecov test coverage](https://codecov.io/gh/jkbest2/spatq/branch/master/graph/badge.svg)](https://codecov.io/gh/jkbest2/spatq?branch=master)
<!-- badges: end -->

This package provides a [Template Model Builder](https://github.com/kaskr/adcomp) program to fit index standardization models that include both fishery-dependent and -independent observations. The goal is to explore the utility of standardizing fishery-dependent observations by allowing for spatially varying catchability.

## Installation

This package has not been released. The most recent version can be installed with:

``` r
devtools::install_github("jkbest2/spatq")
```

## Example

A basic example of generating fake data and fitting a simple model is available in the `tests/testthat/test-spatq.R` file. Currently the `data` and `parameters` arguments to `MakeADFun` must be built by the user. This must be done carefully to avoid mismatched dimensions that can crash your R session.

## Code of conduct

Please note that the 'spatq' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
