estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                          omega = TRUE, epsilon = TRUE,
                          lambda = TRUE, eta = FALSE,
                          phi = FALSE, psi = FALSE)
obj <- make_sim_adfun(repl = 1, sc = "pref",
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
obj_f <- make_sim_adfun(repl = 1, sc = "pref",
                        root_dir = normalizePath(
                          system.file("testdata", package = "spatq")),
                        max_T = 15,
                        spec_estd = estd_f,
                        index_step = 5,
                        ## Can't use Metis if no random effects
                        runSymbolicAnalysis = FALSE,
                        normalize = TRUE,
                        silent = TRUE)
