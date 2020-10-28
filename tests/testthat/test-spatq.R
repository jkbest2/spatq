context("Test TMB model")

test_that("Can get values, gradients, and simulations", {
  fneval <- obj$fn()
  greval <- obj$gr()
  sim <- obj$simulate()

  expect_true(is.finite(fneval))
  expect_true(all(is.finite(greval)))
  expect_true(all(is.finite(unlist(sim))))
})

test_that("Fixed effects model can be fit and sdreported", {
  fneval_f <- obj_f$fn()
  greval_f <- obj_f$gr()
  heeval_f <- obj_f$he()
  fit_f <- fit_spatq(obj_f)
  rep_f <- report_spatq(obj_f)
  sdr_f <- sdreport_spatq(obj_f)

  expect_true(is.finite(fneval_f))
  expect_true(all(is.finite(greval_f)))
  expect_true(all(eigen(heeval_f)$values > 0))

  expect_true(inherits(fit_f, "spatq_fit"))

  expect_true(all(eigen(sdr_f$cov.fixed)$values > 0))
})
