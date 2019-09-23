## Set number of observations, generate observation locations
n_obs <- 1000L
loc <- matrix(runif(2 * n_obs), ncol = 2)
n_years <- 5L
years <- rep(1:n_years, each = 200L)

## Discretize spatial domain into mesh
unit_boundary <- INLA::inla.mesh.segment(matrix(c(0, 0, 1, 1, 0,
                                                  0, 1, 1, 0, 0),
                                                ncol = 2))
mesh <- INLA::inla.mesh.2d(loc,
                           boundary = unit_boundary,
                           offset = c(0.0, 0.25),
                           max.edge = c(0.1, 0.2),
                           min.angle = c(30, 21),
                           cutoff = 0.05)
## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended by
## Finn Lindgren
fem <- INLA::inla.mesh.fem(mesh)
## Rename elements of `fem` so that they are recognized as a `spde` object in
## TMB
names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")

## Make the projection matrix from the vertices to the observation locations
A_spat <- INLA::inla.spde.make.A(mesh, loc)
A_sptemp <- INLA::inla.spde.make.A(mesh, loc, group = years)

## Select observations that are fishery-dependent and generate projection matrix
## to only those locaitons
fishdep <- rbinom(n_obs, 1, 0.75)
A_qspat <- INLA::inla.spde.make.A(mesh, loc, weights = fishdep)
A_qsptemp <- INLA::inla.spde.make.A(mesh, loc, group = years, weights = fishdep)

## Set simulation parameter values. `beta_n` at 0.75 gives roughly 12% zero
## catches, and `beta_w` of 5.0 gives a mean positive catch of about 150. A rho
## of 0.4 means that edges are effectively uncorrelated with the center of the
## unit square.
rho <- 0.4
sigma_spat <- 0.5
kappa <- sqrt(8) / rho
pars_gen <- list(beta_n = rep(0.75, n_years),
                 beta_w = rep(5.0, n_years),
                 lambda_n = rnorm(1),
                 lambda_w = rnorm(1),
                 rho = rho,
                 sigma_spat = sigma_spat,
                 kappa = kappa,
                 tau = sqrt(1 / (4 * pi * kappa ^ 2 * sigma_spat)),
                 sigma_c = 0.25)

## Construct data list that has correct dimensions for use in generating
## simulated data, then do the same with the parameters. The design matrix for
## fixed effects contains only an intercept term for each process.
dat <- list(catch_obs = rep(0, n_obs),
            X_n = model.matrix(~ factor(years) + 0),
            X_w = model.matrix(~ factor(years) + 0),
            A_spat = A_spat,
            A_sptemp = A_sptemp,
            R_n = matrix(rnorm(n_obs), ncol = 1L),
            R_w = matrix(rnorm(n_obs), ncol = 1L),
            A_qspat = A_qspat,
            A_qsptemp = A_qsptemp,
            spde = fem)

## Spatial random effects (`spat_n` and `spat_w`) are set to zero vectors of
## appropriate length, as these are simulated.
pars <- list(beta_n = pars_gen$beta_n,
             beta_w = pars_gen$beta_w,
             omega_n = rep(0.0, mesh$n),
             omega_w = rep(0.0, mesh$n),
             epsilon_n = matrix(0.0, nrow = mesh$n, ncol = n_years),
             epsilon_w = matrix(0.0, nrow = mesh$n, ncol = n_years),
             lambda_n = pars_gen$lambda_n,
             lambda_w = pars_gen$lambda_w,
             phi_n = rep(0.0, mesh$n),
             phi_w = rep(0.0, mesh$n),
             psi_n = matrix(0.0, nrow = mesh$n, ncol = n_years),
             psi_w = matrix(0.0, nrow = mesh$n, ncol = n_years),
             log_kappa = rep(log(pars_gen$kappa), 8),
             log_tau = rep(log(pars_gen$tau), 8),
             log_sigma = log(pars_gen$sigma_c))

map_empty <- list()

## A TMB model object that can be used to simulate
obj <- TMB::MakeADFun(data = dat,
                      parameters = pars,
                      map = map_empty,
                      random = c("omega_n", "omega_w",
                                 "epsilon_n", "epsilon_w",
                                 "phi_n", "phi_w",
                                 "psi_n", "psi_w"),
                      DLL = "spatq",
                      silent = TRUE)

