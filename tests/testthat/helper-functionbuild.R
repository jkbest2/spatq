obj_sim <- make_sim_adfun(repl = 1, sc = "pref",
                          sub_df = data.frame(vessel_idx = 2,
                                              n = 50),
                          root_dir = normalizePath(
                            system.file("testdata", package = "spatq")),
                          max_T = 15,
                          map_pars = c("gamma_n", "gamma_w",
                                       ## "omega_n", "omega_w",
                                       "epsilon_n", "epsilon_w",
                                       "eta_n", "eta_w",
                                       "phi_n", "phi_w",
                                       "psi_n", "psi_w"),
                          runSymbolicAnalysis = TRUE,
                          silent = FALSE)

## fit_sim <- fit_spatq(obj_sim)#, control = list(trace = 2L))
## rep_sim <- report_spatq(obj_sim)
## sdr_sim <- sdreport_spatq(obj_sim)


## par <- obj_sim$env$last.par.best
## ## Use `sapply` because `lapply` doesn't have `USE.NAMES` argument
## par2 <- sapply(unique(names(par)),
##                function(n) unname(par[names(par) == n]),
##                simplify = FALSE, USE.NAMES = TRUE)

## B_true <- read_popstate(1, "simple", system.file("testdata", package = "spatq"))
## I_true <- rescale_index(B_true$pop)

## par(mfrow = c(1, 2))
## plot(par2$beta_n, type = 'l')
## plot(par2$beta_w, type = 'l')

## B_simp <- 25 * exp(par2$beta_n + par2$beta_w)
## I_simp <- rescale_index(B_simp)
## plot(I_simp, type = 'l')
## lines(I_true, col = "red")

## B_est <- 25 * exp(par2$beta_n + mean(par2$omega_n) + par2$beta_w + mean(par2$omega_w))
## I_est <- rescale_index(B_est)
## plot(I_est, type = 'l')
## lines(I_true, col = "red")


## IA_spat <- obj_sim$env$data$IA_spat

## on_ <- IA_spat %*% par2$omega_n
## on <- matrix(on_[1:400], nrow = 20, ncol = 20)
## filled.contour(on, asp = 1)

## ow_ <- IA_spat %*% par2$omega_w
## ow <- matrix(ow_[1:400], nrow = 20, ncol = 20)
## filled.contour(ow, asp = 1)

## pn_ <- IA_spat %*% par2$phi_w
## pn <- matrix(pn_[1:400], nrow = 20, ncol = 20)
## filled.contour(pn, asp = 1)

## pw_ <- IA_spat %*% par2$phi_w
## pw <- matrix(pw_[1:400], nrow = 20, ncol = 20)
## filled.contour(pw, asp = 1)

