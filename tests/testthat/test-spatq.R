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
  ## heeval_f <- obj_f$he()
  fit_f <- fit_spatq(obj_f, NULL, spatq_optcontrol(maxopts = 3))
  rep_f <- report_spatq(obj_f)
  sdr_f <- sdreport_spatq(obj_f)
  ## Want to test that the Hessian is PD *after* getting the parameter
  ## estimates. No expectation that
  ## heeval_f <- obj_f$he(fit_f$par)

  expect_true(is.finite(fneval_f))
  expect_true(all(is.finite(greval_f)))
  ## expect_true(all(eigen(heeval_f)$values > 0))

  expect_true(inherits(fit_f, "spatq_fit"))

  expect_true(all(eigen(sdr_f$cov.fixed)$values > 0))
})

test_that("Fixed effects model with Tweedie observation likelihood works", {
  twf_fn <- twf_obj$fn()
  twf_gr <- twf_obj$gr()
  twf_rep <- report_spatq(twf_obj)
  # tw_sim <- tw_obj$simulate()
  ## tw_fit <- spatq_fit(tw_obj, NULL)

  expect_true(is.finite(twf_fn))
  expect_true(all(is.finite(twf_gr)))
  expect_gt(twf_rep$disp, 0)
  expect_gt(twf_rep$shape, 1)
  expect_lt(twf_rep$shape, 2)
})

test_that("Model with Tweedie observation likelihood works", {
  tw_fn <- tw_obj$fn()
  tw_gr <- tw_obj$gr()
  # tw_sim <- tw_obj$simulate()
  ## tw_fit <- spatq_fit(tw_obj, NULL)

  expect_true(is.finite(tw_fn))
  expect_true(all(is.finite(tw_gr)))
})
