context("Test objective function")

set.seed(394870)

n_obs <- 1000L
loc <- matrix(runif(2 * n_obs), ncol = 2)
mesh <- inla.mesh.2d(loc,
                     offset = c(0.1, 0.2),
                     max.edge = c(0.1, 0.2),
                     min.angle = c(30, 21),
                     cutoff = 0.05)
A <- inla.spde.make.A(mesh, loc)

rho <- 0.4
sigma2 <- 0.5
kappa2 <- 8 / rho^2
tau <- sqrt(1 / (4 * pi * kappa2 * sigma2))
fem <- inla.mesh.fem(mesh)
Q <- tau^2 * (kappa2^2 * fem$c0 + 2 * kappa2 * fem$g1 + fem$g2)
L <- chol(Q)
spat <- backsolve(L, rnorm(mesh$n))
omega <- as.vector(A %*% spat)

sigma_c <- 0.25
## obs <- rlnorm(n_obs, omega, 0.25)
obs <- rlnorm(n_obs, 0, sigma_c)

data <- list(X = matrix(rep(1, n_obs)),
             Y = obs,
             fem = fem,
             A = A)
pars <- list(beta = c(0.0),
             ## log_kappa2 = log(kappa2),
             ## log_tau = log(tau),
             ## spat = spat,
             log_sigma = log(sigma_c))
obj <- TMB::MakeADFun(data = data,
                      parameters = pars,
                      map = list(),
                      ## random = "spat",
                      DLL = "spatq")
obj$fn()
obj$gr()
opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

test_that("Parameters estimated reasonably well", {
  expect_equal(opt$par, unlist(pars),
               tolerance = 1e-1)
})
