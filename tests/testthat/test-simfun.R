context("Test simulation ADFun")

test_that("Can get values, gradients, and simulations", {
  fneval <- obj_fbuild$fn()
  greval <- obj_fbuild$gr()
  sim <- obj_fbuild$simulate()

  expect_true(is.finite(fneval))
  expect_true(all(is.finite(greval)))
  expect_true(all(is.finite(unlist(sim))))
})
