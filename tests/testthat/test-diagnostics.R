context("Test diagnostic functions")

test_that("Null eigenvalue detection works", {
  ## Construct a nearly-singular "covariance" matrix
  n <- 10
  x <- matrix(rnorm(n * n), nrow = n)
  x[, n] <- x[, seq_len(n - 1)] %*% rep(1/9, n - 1) + rnorm(n, 0, 0.001)
  X <- t(x) %*% x

  pars <- rep(0, n)
  names(pars) <- letters[seq_len(n)]
  sdr <- list(par.fixed = pars, cov.fixed = X)

  nulleigs <- fixpar_nulleigs(sdr)
  ## Make sure that the last column is detected as contributing most to singularity.
  expect_equal(nulleigs$sortedpars[length(nulleigs$values), 1], letters[n])
})
