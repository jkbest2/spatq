context("Test utility functions")

test_that("Gathering named vector works", {
  vec <- c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4)
  names(vec) <- c("a", "a", "a", "a", "b", "b", "b", "c", "c", "d")
  glist <- gather_nvec(vec)
  expect_equal(glist,
               list(a = rep(1, 4),
                    b = rep(2, 3),
                    c = rep(3, 2),
                    d = 4))
})
