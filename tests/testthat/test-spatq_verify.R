context("Test input verifications")

# For name checking; missing a component
dat_drop <- dat[-1]
pars_drop <- pars[-3]

test_that("Name checking works", {
  expect_true(verify_spatq_names(dat, pars, map_empty))
  expect_error(verify_spatq_names(dat_drop, pars, map_empty))
  expect_error(verify_spatq_names(dat, pars_drop, map_empty))
})

## For dim checks
dat_bad <- dat
dat_bad$X_n <- dat_bad$X_n[, -2]

pars_bad <- pars
pars_bad$beta_n <- pars_bad$beta_n[-2]

## For dim and map checking
map_epsilon_n_bad <- list(epsilon_n = matrix(NA,
                                             nrow(pars$epsilon_n),
                                             ncol(pars$epsilon_n)))
map_epsilon_n <- map_epsilon_n_bad
map_epsilon_n$log_kappa <- factor(c(1, 2, NA, 3, 4, 5, 6, 7))
map_epsilon_n$log_tau <- factor(c(1, 2, NA, 3, 4, 5, 6, 7))

map_xi_bad <- list(gamma_n = factor(rep(NA, length(pars$gamma_n))))
map_xi_bad2 <- map_xi_bad
map_xi_bad2$log_xi <- factor(c(1, 2, 3, 4))
map_xi <- map_xi_bad2
map_xi$log_xi[1] <- NA

test_that("Dimension checking works", {
  expect_true(verify_spatq_dims(dat, pars, map_empty))
  expect_true(verify_spatq_dims(dat, pars, map_epsilon_n))
  expect_error(verify_spatq_dims(dat, pars_bad, map_empty))
  expect_error(verify_spatq_dims(dat_bad, pars, map_empty))
})

test_that("Map checking works", {
  expect_true(verify_spatq_map(pars, map_empty))
  expect_true(verify_spatq_map(pars, map_epsilon_n))
  expect_error(verify_spatq_map(pars, map_epsilon_n_bad))
  expect_true(verify_spatq_map(pars, map_xi, FALSE))
  expect_error(verify_spatq_map(pars, map_xi_bad))
  expect_error(verify_spatq_map(pars, map_xi_bad2))
})

