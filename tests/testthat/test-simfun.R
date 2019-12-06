context("Test simulation ADFun")

test_that("Can get values, gradients, and simulations", {
  fneval <- obj_sim$fn()
  greval <- obj_sim$gr()
  sim <- obj_sim$simulate()

  expect_true(is.finite(fneval))
  expect_true(all(is.finite(greval)))
  expect_true(all(is.finite(unlist(sim))))
})
