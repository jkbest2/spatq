obj_fbuild <- make_sim_adfun(repl = 1, sc = "simple",
                             sub_df = data.frame(vessel_idx = 2,
                                                 n = 1000),
                             root_dir = normalizePath(
                               system.file("testdata", package = "spatq")),
                             max_T = 5)

## fit <- fit_spatq(obj_fbuild)
## rep <- report_spatq(obj_fbuild)
## sdr <- sdreport_spatq(obj_fbuild)
