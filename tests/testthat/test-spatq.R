context("Test TMB model")

## Set number of observations, generate observation locations
n_obs <- 1000L
loc <- matrix(runif(2 * n_obs), ncol = 2)

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
A <- INLA::inla.spde.make.A(mesh, loc)

## Set simulation parameter values. `beta_n` at 0.75 gives roughly 12% zero
## catches, and `beta_w` of 5.0 gives a mean positive catch of about 150. A rho
## of 0.4 means that edges are effectively uncorrelated with the center of the
## unit square.
rho <- 0.4
sigma_spat <- 0.5
kappa <- sqrt(8) / rho
pars_gen <- list(beta_n = 0.75,
                 beta_w = 5.0,
                 rho = rho,
                 sigma_spat = sigma_spat,
                 kappa = kappa,
                 tau = sqrt(1 / (4 * pi * kappa ^ 2 * sigma_spat)),
                 sigma_c = 0.25)

## Construct data list that has correct dimensions for use in generating
## simulated data, then do the same with the parameters. The design matrix for
## fixed effects contains only an intercept term for each process.
data <- list(catch_obs = rep(0, n_obs),
             X_n = matrix(rep(1, n_obs)),
             X_w = matrix(rep(1, n_obs)),
             spde = fem,
             A = A)
## Spatial random effects (`spat_n` and `spat_w`) are set to zero vectors of
## appropriate length, as these are simulated.
pars <- list(beta_n = pars_gen$beta_n,
             beta_w = pars_gen$beta_w,
             log_kappa = rep(log(pars_gen$kappa), 2),
             log_tau = rep(log(pars_gen$tau), 2),
             log_sigma = log(pars_gen$sigma_c),
             spat_n = rep(0.0, mesh$n),
             spat_w = rep(0.0, mesh$n))

## A TMB model object that can be used to simulate
obj <- TMB::MakeADFun(data = data,
                      parameters = pars,
                      map = list(),
                      random = c("spat_n", "spat_w"),
                      DLL = "spatq")

test_that("Simulation study - long form", {
  ## Don't try to run 6 hours of tests on Travis, it's probably frowned upon
  skip_on_travis()

  ## 50 seems like a reasonable number of simulation, currently takes around 6
  ## hours to complete. Use a global variable to keep track of progress for the
  ## progress bar. Seems hacky, but not too bad.
  n_repl <- 5L
  rep_cnt <- 0L
  ## pb <- txtProgressBar(min = 0, max = n_repl, style = 3)
  message("Started: ", Sys.time())
  sim <- replicate(n_repl, {
    simdata <- obj$simulate(par = unlist(pars), complete = TRUE)
    obj2 <- TMB::MakeADFun(simdata, pars,
                           random = c("spat_n", "spat_w"),
                           DLL = "spatq", silent = TRUE)
    start_time <- Sys.time()
    fit <- tryCatch({
      opt <- optim(obj2$par, obj2$fn, obj2$gr,
                   method = "BFGS",
                   control = list(maxit = 100L))
      ## Keep the entire vector of parameters (including spatial random
      ## effects). Also return the simulated data for comparison, and the
      ## optimization result to check for e.g. convergence and number of
      ## function evaluations.
      list(simdata = simdata,
           pars = obj2$env$last.par.best,
           opt = opt)
    }, error = function(cnd) {
      ## Return the same list, but include NAs where appropriate. Having the
      ## simulations that fail might be helpful for identifying where and why
      ## they fail.
      val <- list(simdata = simdata,
                  pars = NA,
                  opt = NA)
      attr(val, "error") <- as.character(cnd)
      val
    })
    fit$elapsed_time <- Sys.time() - start_time
    rep_cnt <<- rep_cnt + 1L
    ## setTxtProgressBar(pb, rep_cnt)
    message("Finished ", rep_cnt, ": ", Sys.time())
    fit
  }, simplify = FALSE) # Always return a list

  ## Check that not too many fail to converge. Counting failure here as both
  ## erroring out and non-0 convergence code. 95% success sounds like a good
  ## starting goal.
  fail_idx <- vapply(sim,
                     function(s) {
                       is.na(s$opt) || s$opt$convergence > 0L
                     }, FUN.VALUE = TRUE)
  expect_lt(mean(fail_idx), 0.05)

  ## Test that generative values are within some quantile range of the fitted
  ## values, using only the models that (are claimed to have) converged.
  sim_ests <- Reduce(rbind, lapply(sim[!fail_idx], function(s) s$pars[1:7]))
  sim_q10 <- apply(sim_est, 2L, quantile, prob = 0.1)
  sim_q90 <- apply(sim_est, 2L, quantile, prob = 0.9)
  expect_true(all(sim_q10 < obj$par))
  expect_true(all(obj$par < sim_q90))

  ## Test that none of the converged models include parameter values that are
  ## exactly equal to their initial values, indicating the parameter never
  ## moved.
  expect_false(any(obj2$par == t(sim_est)))

  ## Broad performance metric - at least warn me if fits start taking longer
  sim_eltime <- vapply(sim[!fail_idx],
                       function(s) as.double(s$elapsed_time, units = "mins"), 1.0)
  ## Somewhat arbitrary, but maybe worry if the median successful fit is over 15
  ## minutes and the 90th percentile is over 20 minutes?
  expect_lt(median(sim_eltime), 15.0)
  expect_lt(quantile(sim_eltime, prob = 0.9), 20)
}
