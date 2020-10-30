context("Test setup functions")

test_that("Single parameter update works", {
  mapvec1 <- c(1, 2, NA, NA, 1, 2, NA, NA)
  mapvec2 <- c(1, 2, NA, NA, 3, 4, NA, NA)
  mapvec3 <- c(1, 2, NA, NA, 3, 4, NA, NA, 5, 6)
  newvals1 <- c(100, 200)
  newvals2 <- seq(100, 800, by = 100)
  currvals1 <- 1:8
  currvals2 <- 1:10

  ## Straightforward cases
  expect_equal(update_onepar(currvals1),
               currvals1)
  expect_equal(update_onepar(currvals1, newvals2),
               newvals2)
  expect_equal(update_onepar(currvals1, newvals1),
               c(100, 200, 3, 4, 5, 6, 7, 8))
  expect_equal(update_onepar(currvals1, newvals1, mapvec1),
               c(100, 200, 3, 4, 100, 200, 7, 8))
  expect_equal(update_onepar(currvals1, newvals1, mapvec2),
               c(100, 200, 3, 4, 5, 6, 7, 8))
  expect_equal(update_onepar(currvals2, newvals1, mapvec3),
               c(100, 200, 3, 4, 5, 6, 7, 8, 9, 10))
  expect_error(update_onepar(currvals1, newvals1, mapvec3),
               "Length of currvals must be greater than mapvec")


})
