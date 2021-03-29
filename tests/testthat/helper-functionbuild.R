estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                          omega = TRUE, epsilon = TRUE,
                          lambda = TRUE, eta = FALSE,
                          phi = FALSE, psi = FALSE)
obj <- make_sim_adfun(study = "qdevscaling",
                      repl = 1,
                      opmod = 1,
                      sub_df = data.frame(vessel_idx = 2,
                                          n = 500),
                      root_dir = normalizePath(
                        system.file("testdata", package = "spatq")),
                      max_T = 15,
                      spec_estd = estd,
                      index_step = 5,
                      runSymbolicAnalysis = TRUE,
                      normalize = TRUE,
                      silent = TRUE)

## Fixed-effect only model
estd_f <- specify_estimated()
obj_f <- make_sim_adfun(study = "qdevscaling",
                        repl = 1,
                        opmod = 1,
                        root_dir = normalizePath(
                          system.file("testdata", package = "spatq")),
                        max_T = 15,
                        spec_estd = estd_f,
                        index_step = 5,
                        ## Can't use Metis if no random effects
                        runSymbolicAnalysis = FALSE,
                        normalize = TRUE,
                        silent = TRUE)

## Tweedie with fixed effects
twf_estd <- specify_estimated(beta = TRUE, lambda = TRUE,
                              obs_lik = 1L)
twf_setup <- spatq_simsetup(study = "qdevscaling",
                            repl = 1,
                            opmod = 1,
                            root_dir = normalizePath(
                              system.file("testdata", package = "spatq")),
                            max_T = 5,
                            index_step = 50,
                            spec_estd = twf_estd)
twf_obj <- spatq_obj(twf_setup)
twf_fit <- fit_spatq(twf_obj)

## Tweedie with random effects
tw_estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                             omega = TRUE, epsilon = TRUE,
                             lambda = TRUE, eta = FALSE,
                             phi = TRUE, psi = FALSE,
                             kappa_map = c(1, NA, 1, NA, 1, NA, NA, NA),
                             obs_lik = 1L)
tw_setup <- spatq_simsetup(study = "qdevscaling",
                           repl = 1,
                           opmod = 1,
                           root_dir = normalizePath(
                             system.file("testdata", package = "spatq")),
                           max_T = 15,
                           index_step = 5,
                           spec_estd = tw_estd,
                           init_fixef = TRUE)
tw_obj <- spatq_obj(tw_setup)
