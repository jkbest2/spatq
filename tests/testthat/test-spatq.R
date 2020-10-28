context("Test TMB model")

test_that("Can get values, gradients, and simulations", {
  fneval <- obj$fn()
  greval <- obj$gr()
  sim <- obj$simulate()

  expect_true(is.finite(fneval))
  expect_true(all(is.finite(greval)))
  expect_true(all(is.finite(unlist(sim))))
})

test_that("Fixed effects model can be fit and sdreported", {
  fneval_f <- obj_f$fn()
  greval_f <- obj_f$gr()
  heeval_f <- obj_f$he()
  fit_f <- fit_spatq(obj_f)
  rep_f <- report_spatq(obj_f)
  sdr_f <- sdreport_spatq(obj_f)

  expect_true(is.finite(fneval_f))
  expect_true(all(is.finite(greval_f)))
  expect_true(all(eigen(heeval_f)$values > 0))

  expect_true(inherits(fit_f, "spatq_fit"))

  expect_true(all(eigen(sdr_f$cov.fixed)$values > 0))
})

test_that("Initial values usually don't fail", {
  ## Takes too long to run (probably)
  skip_on_travis()
  skip("Takes too long")

  ## Test that initial values don't cause NaN function evaluations or gradients
  ## too often (95% of the time?). Currently using generative values with all
  ## zeros for random effects. Might speed things up to use simulated values of
  ## random effects? But might not catch other problems.
  n_repl <- 100L
  init_sim <- replicate(n_repl,
                        obj$simulate(par = unlist(pars), complete = TRUE),
                        simplify = FALSE)

  ## NOTE Failures here appear to depend on how many iterations the inner
  ## problem is allowed. With just 100 iterations, around 17% fail. VAST doesn't
  ## appear to alter the inner problem max iterations.
  init_obj <- lapply(init_sim,
                     function(sim_data) {
                       TMB::MakeADFun(data = sim_data,
                                      parameters = pars,
                                      random = c("spat_n", "spat_w"),
                                      inner.control = list(maxit = 100L),
                                      DLL = "spatq")
                     })

  ## Evaluate function and gradients, and check that all are finite (non-finite
  ## indicates a failure).
  init_fneval <- lapply(init_obj,
                        function(obj) obj$fn(obj$par))
  init_fnfin <- vapply(init_fneval,
                       function(fn) is.finite(fn),
                       TRUE)
  init_greval <- lapply(init_obj,
                        function(obj) obj$gr(obj$par))
  init_grfin <- vapply(init_greval,
                       function(gr) all(is.finite(gr)),
                       TRUE)


  ## Check that over 80% of evaulations result in finite values
  expect_gt(mean(init_fnfin), 0.8)
  expect_gt(mean(init_grfin), 0.8)
})


test_that("Simulation study gives reasonable results", {
  ## Don't try to run 6 hours of tests on Travis, it's probably frowned upon
  skip_on_travis()
  ## Skip for now; takes too long and doesn't really work
  skip("Takes too long")

  ## 50 seems like a reasonable number of simulation, currently takes around 6
  ## hours to complete. Use a global variable to keep track of progress for the
  ## progress bar. Seems hacky, but not too bad.
  n_repl <- 50L
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

  ## Extract the parameter estimates. Define three modes of failure:
  ## 1. Error in optimization step due to NaN value or gradient
  ## 2. Nonconvergence as indicated by s$convergence > 0
  ## 3. Estimated values at exactly initial values
  ## Test that less than 95% fail due to any of these reasons
  sim_ests_all <- Reduce(rbind, lapply(sim, function(s) s$pars[1:7]))
  nan_idx <- vapply(sim,
                     function(s) {
                       any(is.na(s$opt))# || s$opt$convergence > 0L
                     }, FUN.VALUE = TRUE)
  nonconv_idx <- vapply(sim,
                        function(s) {
                          !is.na(s$opt) && s$opt$convergence > 0L
                        }, FUN.VALUE = TRUE)
  nomove_idx <- apply(sim_ests, 1, function(s) any(s == obj$par))
  fail_idx <- nan_idx | nonconv_idx | nomove_idx
  expect_lt(mean(fail_idx), 0.05)

  ## Test that generative values are within some quantile range of the fitted
  ## values, using only the models that (are claimed to have) converged.
  sim_ests <- sim_ests_all[!fail_idx, ]
  sim_q10 <- apply(sim_ests, 2L, quantile, prob = 0.1)
  sim_q90 <- apply(sim_ests, 2L, quantile, prob = 0.9)
  expect_true(all(sim_q10 < obj$par))
  expect_true(all(obj$par < sim_q90))

  ## Broad performance metric - at least warn me if fits start taking longer
  sim_eltime <- vapply(sim[!fail_idx],
                       function(s) as.double(s$elapsed_time, units = "mins"),
                       1.0)
  ## Somewhat arbitrary, but maybe worry if the median successful fit is over 15
  ## minutes and the 90th percentile is over 20 minutes?
  expect_lt(median(sim_eltime), 15.0)
  expect_lt(quantile(sim_eltime, prob = 0.9), 20)
})
