obj_fbuild <- make_sim_adfun(repl = 1, sc = "simple",
                             sub_df = data.frame(vessel_idx = c(1, 2),
                                                 n = c(100, 100)),
                             root_dir = normalizePath(
                               system.file("testdata", package = "spatq")))

