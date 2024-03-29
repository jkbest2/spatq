### Construct all pieces of a spatq model manually to compare with the generated
### pieces
## Set number of observations, generate observation locations
n_obs_year <- 1500L
n_years <- 15L
n_obs <- n_obs_year * n_years
loc <- matrix(runif(2 * n_obs), ncol = 2)
years <- rep(1:n_years, each = n_obs_year)

## Generate index locations
## TODO Extract index grid-making to a function
index_step <- 0.05
index_start <- index_step / 2
index_end <- 1 - index_start
index_vec <- seq(index_start, index_end, index_step)
index_coords <- expand.grid(
  s1 = index_vec,
  s2 = index_vec,
  year = 1:n_years
)
n_index <- nrow(index_coords)
index_years <- index_coords$year
index_loc <- as.matrix(index_coords[c("s1", "s2")])

Ih <- rep_len(index_step^2, n_index)

## Discretize spatial domain into mesh
unit_boundary <- INLA::inla.mesh.segment(matrix(c(
  0, 0, 1, 1, 0,
  0, 1, 1, 0, 0
),
ncol = 2
))
mesh <- INLA::inla.mesh.2d(loc,
  boundary = unit_boundary,
  offset = c(0.0, 0.25),
  max.edge = c(0.1, 0.2),
  min.angle = c(30, 21),
  cutoff = 0.05
)
## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended by
## Finn Lindgren
## fem <- INLA::inla.mesh.fem(mesh)
## Rename elements of `fem` so that they are recognized as a `spde` object in
## TMB
## names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")
fem <- generate_aniso_fem(mesh)

## Make the projection matrix from the vertices to the observation locations
A_spat <- INLA::inla.spde.make.A(mesh, loc)
A_sptemp <- INLA::inla.spde.make.A(mesh, loc, group = years)
## And the index projection matrices
IA_spat <- INLA::inla.spde.make.A(mesh, index_loc)
IA_sptemp <- INLA::inla.spde.make.A(mesh, index_loc, group = index_years)

## Select observations that are fishery-dependent and generate projection matrix
## to only those locations
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
pars_gen <- list(
  beta_n = rep(0.75, n_years),
  beta_w = rep(5.0, n_years),
  gamma_n = rnorm(10),
  gamma_w = rnorm(10),
  lambda_n = rnorm(1),
  lambda_w = rnorm(1),
  eta_n = rnorm(10),
  eta_w = rnorm(10),
  rho = rho,
  sigma_spat = sigma_spat,
  kappa = kappa,
  tau = sqrt(1 / (4 * pi * kappa^2 * sigma_spat)),
  H_pars = c(log(0), 0),
  sigma_c = 0.25
)

## Construct data list that has correct dimensions for use in generating
## simulated data, then do the same with the parameters. The design matrix for
## fixed effects contains only an intercept term for each process.
dat <- list(
  catch_obs = rep(0, n_obs),
  area_swept = rep(1, n_obs),
  X_n = model.matrix(~ factor(years) + 0),
  X_w = model.matrix(~ factor(years) + 0),
  IX_n = model.matrix(~ factor(index_years) + 0),
  IX_w = model.matrix(~ factor(index_years) + 0),
  Z_n = model.matrix(~ factor(sample(1:10, n_obs, TRUE)) + 0),
  Z_w = model.matrix(~ factor(sample(1:10, n_obs, TRUE)) + 0),
  IZ_n = model.matrix(~ factor(sample(1:10, n_index, TRUE)) + 0),
  IZ_w = model.matrix(~ factor(sample(1:10, n_index, TRUE)) + 0),
  A_spat = A_spat,
  A_sptemp = A_sptemp,
  IA_spat = IA_spat,
  IA_sptemp = IA_sptemp,
  Ih = Ih,
  R_n = matrix(rnorm(n_obs), ncol = 1L),
  R_w = matrix(rnorm(n_obs), ncol = 1L),
  V_n = model.matrix(~ factor(sample(1:10, n_obs, TRUE)) + 0),
  V_w = model.matrix(~ factor(sample(1:10, n_obs, TRUE)) + 0),
  A_qspat = A_qspat,
  A_qsptemp = A_qsptemp,
  spde = fem,
  proc_switch = rep(TRUE, 6),
  obs_lik = 0L,
  norm_flag = FALSE,
  incl_data = TRUE
)

## Spatial random effects (`spat_n` and `spat_w`) are set to zero vectors of
## appropriate length, as these are simulated.
pars <- list(
  beta_n = pars_gen$beta_n,
  beta_w = pars_gen$beta_w,
  gamma_n = pars_gen$gamma_n,
  gamma_w = pars_gen$gamma_w,
  omega_n = rep(0.0, mesh$n),
  omega_w = rep(0.0, mesh$n),
  epsilon_n = matrix(0.0, nrow = mesh$n, ncol = n_years),
  epsilon_w = matrix(0.0, nrow = mesh$n, ncol = n_years),
  lambda_n = pars_gen$lambda_n,
  lambda_w = pars_gen$lambda_w,
  eta_n = pars_gen$eta_n,
  eta_w = pars_gen$eta_w,
  phi_n = rep(0.0, mesh$n),
  phi_w = rep(0.0, mesh$n),
  psi_n = matrix(0.0, nrow = mesh$n, ncol = n_years),
  psi_w = matrix(0.0, nrow = mesh$n, ncol = n_years),
  log_xi = rep(0.0, 4L),
  log_kappa = rep(log(pars_gen$kappa), 8),
  log_tau = rep(log(pars_gen$tau), 8),
  H_pars = c(log(0), 0),
  obs_lik_pars = log(pars_gen$sigma_c)
)
attr(pars, "map_lambda") <- FALSE

map_empty <- list()

verify_spatq(dat, pars, map_empty)
