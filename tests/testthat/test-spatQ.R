context("Test objective function")

test_that("Compiled model runs and converges", {

  set.seed(394870)

  n_obs <- 10000L
  loc <- matrix(runif(2 * n_obs), ncol = 2)
  mesh <- INLA::inla.mesh.2d(loc,
                       offset = c(0.1, 0.2),
                       max.edge = c(0.1, 0.2),
                       min.angle = c(30, 21),
                       cutoff = 0.05)
  A <- INLA::inla.spde.make.A(mesh, loc)

  rho <- 0.4
  sigma2 <- 0.5
  kappa2 <- 8 / rho^2
  tau <- sqrt(1 / (4 * pi * kappa2 * sigma2))
  fem <- INLA::inla.mesh.fem(mesh)
  Q <- tau^2 * (kappa2^2 * fem$c0 + 2 * kappa2 * fem$g1 + fem$g2)
  L <- chol(Q)
  spat <- backsolve(L, rnorm(mesh$n))
  omega <- as.vector(A %*% spat)

  ## Change names to M0, M1, M2
  names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")

  sigma_c <- 0.25
  ## obs <- rlnorm(n_obs, omega, 0.25)
  obs <- rlnorm(n_obs, 0, sigma_c)
  obs[sample(length(obs), 2000)] <- 0

  data <- list(catch_obs = obs,
               X_n = matrix(rep(1, n_obs)),
               X_w = matrix(rep(1, n_obs)),
               spde = fem,
               A = A)
  pars <- list(beta_n = c(0.0),
               beta_w = c(0.0),
               log_kappa = rep(0.5 * log(kappa2), 2),
               log_tau = rep(log(tau), 2),
               log_sigma = log(sigma_c),
               spat_n = rep(0.0, mesh$n),
               spat_w = rep(0.0, mesh$n))
  obj <- TMB::MakeADFun(data = data,
                        parameters = pars,
                        map = list(),
                        random = c("spat_n", "spat_w"),
                        DLL = "spatq")
  opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

  expect_equal(opt$convergence, 0L)
})
