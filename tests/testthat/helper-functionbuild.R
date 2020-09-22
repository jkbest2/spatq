estd <- specify_estimated(beta = TRUE,
                          gamma = FALSE,
                          omega = TRUE,
                          epsilon = TRUE,
                          lambda = TRUE,
                          eta = FALSE,
                          phi = FALSE,
                          psi = FALSE)
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

## fit_sim <- fit_spatq(obj)#, control = list(trace = 2L))
## rep_sim <- report_spatq(obj)
## sdr_sim <- sdreport_spatq(obj)
