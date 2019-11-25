obj_fbuild <- make_sim_adfun(repl = 1, sc = "simple",
                             sub_df = data.frame(vessel_idx = c(1, 2),
                                                 n = c(400, 500)),
                             root_dir = normalizePath(
                               system.file("testdata", package = "spatq")),
                             max_T = 15,
                             runSymbolicAnalysis = FALSE)

## fit <- fit_spatq(obj_fbuild, control = list(trace = 2L))
## rep <- report_spatq(obj_fbuild)
## sdr <- sdreport_spatq(obj_fbuild)
