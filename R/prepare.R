read_catch <- function(repl, sc) {
  sc %in% c("naive", "simple", "scaled", "shared")
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <-str_pad(repl, 2, pad = "0")
  flnm <- paste0("repl_", repl_str, "/catch_", repl_str, "_", sc, ".csv")
 read_csv(flnm,
                  col_types =cols(
                    time =col_integer(),
                    vessel_idx =col_integer(),
                    loc_idx =col_integer(),
                    coordinates =col_character(),
                    effort =col_double(),
                    catch_biomass =col_double())) %>%
   left_join(get_coordref(), by = "loc_idx")
}

subsample_catch <- function(catch_df, n_df = NULL) {
  if (is.null(n_df)) return(catch_df)
  nms <- names(n_df)
  "n" %in% names(n_df) || stop("Need column \`n\`")
  nms_join <- nms[nms != "n"]
  catch_df %>%
    nest(data = c(-vessel_idx, -time)) %>%
   left_join(n_df, by = nms_join) %>%
   mutate(n = coalesce(n, map_dbl(data, nrow))) %>%
   mutate(data = map2(data, n, sample_n)) %>%
    select(-n) %>%
    unnest(data)
}

read_popstate <- function(repl, sc) {
  sc %in% c("naive", "simple", "scaled", "shared")
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <-str_pad(repl, 2, pad = "0")
  flnm <- paste0("repl_", repl_str, "/popstate_", repl_str, "_", sc, ".csv")
  pop <-read_csv(flnm, col_types =cols(pop =col_double()))
 mutate(pop, time = seq_along(pop))
}

domain_boundary <- function() {
  boundary <- INLA::inla.mesh.segment(matrix(c(0, 0,   100, 100, 0,
                                               0, 100, 100, 0,   0),
                                             ncol = 2))
}

loc_grid <- function(step = 1.0) {
  as.matrix(expand.grid(s1 = seq(0.5, 99.5, step),
                        s2 = seq(0.5, 99.5, step)))
}

generate_mesh <- function() {
  ## Discretize spatial domain into mesh
  boundary <- domain_boundary()
  loc <- loc_grid(2.0)
  INLA::inla.mesh.2d(loc,
                     boundary = boundary,
                     offset = c(5.0, 25.0),
                     max.edge = c(3 * sqrt(2), 6 * sqrt(2)),
                     min.angle = c(30, 21),
                     cutoff = 0.5)
}

generate_fem <- function(mesh) {
  ## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended
  ## by Finn Lindgren
  fem <- INLA::inla.mesh.fem(mesh)
  ## Rename elements of `fem` so that they are recognized as a `spde` object in
  ## TMB
  names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")
  fem
}

generate_projection <- function(mesh,
                                data_df,
                                vessel_idx = NULL,
                                group = NULL) {
  if (is.null(vessel_idx)) {
    wts <- NULL
  } else {
    wts <- data_df$vessel_idx %in% vessel_idx
  }
  loc <- as.matrix(data_df[, c("s1", "s2")])
  A <- INLA::inla.spde.make.A(mesh, loc, group = group, weights = wts)
  Matrix::drop0(A)
}

generate_empty_projection <- function(mesh, data_df, group = 1, ...) {
  ## Matrix::Matrix(0,
  ##                nrow = nrow(data_df),
  ##                ncol = mesh$n * length(unique(group)))
  INLA::inla.spde.make.A(mesh, group = 1:25)
}

parse_coords <- function(coord_tuple) {
  coord_str <- str_sub(coord_tuple, start = 2L, end = -2L)
  coords <- as.numeric(str_split_fixed(coord_str, ", ", 2))
  coords
}

get_coordref <- function() {
  if (!file.exists("coordref.csv")) {
    stop("The file \'coordref.csv\' is not in the current working directory")
  }
 read_csv("coordref.csv",
                  col_types =cols(
                    loc_idx =col_integer(),
                    s1 =col_double(),
                    s2 =col_double()))
}

create_index_df <- function(step = 1.0, T = 25) {
  loc_df <- as_tibble(loc_grid(step)) %>%
   mutate(loc_idx = seq_along(s1))
  index_df <- cross_df(list(loc_idx = loc_df$loc_idx, time = seq_len(T))) %>%
    full_join(loc_df, by = "loc_idx") %>%
    select(-loc_idx)
  attr(index_df, "step") <- step
  index_df
}

prepare_data <- function(catch_df, index_df, mesh, fem) {
  T <- length(unique(catch_df$time))

  dat <- list(catch_obs = catch_df$catch_biomass,
              area_swept = catch_df$effort,
              X_n = model.matrix(~ factor(time) + 0, data = catch_df),
              X_w = model.matrix(~ factor(time) + 0, data = catch_df),
              IX_n = model.matrix(~ factor(time) + 0, data = index_df),
              IX_w = model.matrix(~ factor(time) + 0, data = index_df),
              Z_n = matrix(0, nrow = nrow(catch_df), ncol = 1),
              Z_w = matrix(0, nrow = nrow(catch_df), ncol = 1),
              IZ_n = matrix(0, nrow = nrow(index_df), ncol = 1),
              IZ_w = matrix(0, nrow = nrow(index_df), ncol = 1),
              A_spat = generate_projection(mesh, catch_df),
              A_sptemp = generate_projection(mesh, catch_df,
                                           group = catch_df$time),
              ## A_sptemp = generate_empty_projection(mesh, catch_df,
              ##                                      group = catch_df$time),
              IA_spat = generate_projection(mesh, index_df),
              IA_sptemp = generate_projection(mesh, index_df,
                                              group = index_df$time),
              ## IA_sptemp = generate_empty_projection(mesh, index_df,
              ##                                       group = catch_df$time),
              Ih = rep_len(attr(index_df, "step")^2, nrow(index_df)),
              R_n = model.matrix(~ factor(vessel_idx), data = catch_df),
              R_w = model.matrix(~ factor(vessel_idx), data = catch_df),
              V_n = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              V_w = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              A_qspat = generate_projection(mesh, catch_df, vessel_idx = 2),
              A_qsptemp = generate_projection(mesh, catch_df,
                                            vessel_idx = 2,
                                            group = catch_df$time),
              ## A_qsptemp = generate_empty_projection(mesh, catch_df,
              ##                                       group = catch_df$time),
              spde = fem)
  ## Store the number of years as an attribute
  attr(dat, "T") <- T
  dat
}

pars_data <- function(data_el, T = 1) {
  if (inherits(data_el, "inla.mesh")) {
    n <- data_el$n
    if (T > 1) {
      par_el <- matrix(0.0, nrow = n, ncol = T)
    } else {
      par_el <- rep_len(0.0, n)
    }
  } else if (is.null(dim(data_el))) {
    par_el <- 0.0
  } else {
    par_el <- rep_len(0.0, ncol(data_el))
  }
  par_el
}

pars_kappa <- function(rho) {
  sqrt(8) / rho
}

pars_tau <- function(sig2, rho) {
  den <- sig2 * 4 * pi * pars_kappa(rho)^2
  sqrt(1 / den)
}

prepare_pars <- function(data, mesh) {
  T <- attr(data, "T")
  pars <- list(beta_n = pars_data(data$X_n),
               beta_w = pars_data(data$X_w),
               gamma_n = pars_data(data$Z_n),
               gamma_w = pars_data(data$Z_w),
               omega_n = pars_data(mesh),
               omega_w = pars_data(mesh),
               epsilon_n = pars_data(mesh, T),
               epsilon_w = pars_data(mesh, T),
               lambda_n = pars_data(data$R_n),
               lambda_w = pars_data(data$R_w),
               eta_n = pars_data(data$V_n),
               eta_w = pars_data(data$V_w),
               phi_n = pars_data(mesh),
               phi_w = pars_data(mesh),
               psi_n = pars_data(mesh, T),
               psi_w = pars_data(mesh, T),
               log_xi = rep(0.0, 4L),
               log_kappa = rep(log(pars_kappa(50)), 8),
               log_tau = rep(log(pars_tau(1.0, 50)), 8),
               log_sigma = log(1.0))
}

spat_par_idx <- function(par_name) {
  spat_pars <- c("omega_n", "omega_w",
                 "epsilon_n", "epsilon_w",
                 "phi_n", "phi_w",
                 "psi_n", "psi_w")
  which(par_name == spat_pars)
}

re_par_idx <- function(par_name) {
  re_pars <- c("gamma_n", "gamma_w",
               "eta_n", "eta_w")
  which(par_name == re_pars)
}

prepare_map <- function(pars, map_list) {
  map <- lapply(pars, function(par) par[] <- factor(seq_along(par)))
  for (par in map_list) {
    map[[par]][] <- NA
    sp_idx <- spat_par_idx(par)
    if (length(sp_idx) > 0) {
      map$log_kappa[sp_idx] <- NA
      map$log_tau[sp_idx] <- NA
    }
    re_idx <- re_par_idx(par)
    if (length(re_idx) > 0) {
      map$log_xi[re_idx] <- NA
    }
  }
  map
}

prepare_random <- function(map) {
  re_pars <- c("gamma_n", "gamma_w",
               "omega_n", "omega_w",
               "epsilon_n", "epsilon_w",
               "eta_n", "eta_w",
               "phi_n", "phi_w",
               "psi_n", "psi_w")
  is_mapd <- unlist(lapply(re_pars, function(par) all(is.na(map[[par]]))))
  re_pars[!is_mapd]
}

prepare_adfun <- function(data, parameters, map, random,
                          ..., silent = TRUE) {
  verify_spatq(data, parameters, map)
  TMB::MakeADFun(data = data,
                 parameters = parameters,
                 map = map,
                 random = random,
                 DLL = "spatq",
                 ...,
                 silent = silent)
}

make_sim_adfun <- function(repl, sc, sub_df = NULL) {
  ## Read in data
  catch_df <- read_catch(repl, sc)
  ## Subset observations
  catch_df <- subsample_catch(catch_df, sub_df)

  ## Create index integration reference
  index_df <- create_index_df(step = 1, T = 25)

  ## Discretize space
  mesh <- generate_mesh()
  fem <- generate_fem(mesh)

  data <- prepare_data(catch_df, index_df, mesh, fem)
  parameters <- prepare_pars(data, mesh)
  map <- prepare_map(parameters,
                     map_list = c("gamma_n", "gamma_w",
                                  "epsilon_n", "epsilon_w",
                                  "eta_n", "eta_w",
                                  "psi_n", "psi_w"))
  random <- prepare_random(map)
  print(nrow(data$IX_n))
  print(nrow(data$IX_w))
  print(nrow(data$IA_spat))
  print(nrow(data$IA_sptemp))
  stopifnot(expr = {
    nrow(data$IX_n) == 250000
    nrow(data$IX_w) == 250000
    nrow(data$IA_spat) == 250000
    nrow(data$IA_sptemp) == 250000
  })
  prepare_adfun(data, parameters, map, random)
}
