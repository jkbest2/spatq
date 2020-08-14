##' Read in a replicate of data simulated by FisherySim.jl
##'
##' Requires that catch files are e.g. \code{repl_01/catch_simple.csv}.
##' @title Read simulated data
##' @param repl Replicate number
##' @param sc Scenario; one of "pref", "spat", or "combo"
##' @param root_dir Directory to work in; defaults to current working directory
##' @return A \code{tibble} with catch observations
##' @author John Best
##' @importFrom dplyr %>%
##' @export
read_catch <- function(repl, sc, root_dir = ".") {
  check_scenario(sc)
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <- stringr::str_pad(repl, 2, pad = "0")
  repl_dir <- paste0("repl_", repl_str)
  sc_file <- paste0("catch_", repl_str, "_", sc, ".csv")
  flnm <- file.path(root_dir, repl_dir, sc_file)
  catch_df <- readr::read_csv(flnm,
                              col_types = readr::cols(
                                time = readr::col_integer(),
                                vessel_idx = readr::col_integer(),
                                loc_idx = readr::col_integer(),
                                coordinates = readr::col_character(),
                                effort = readr::col_double(),
                                catch_biomass = readr::col_double())) %>%
    dplyr::left_join(get_coordref(root_dir), by = "loc_idx")
  ## Number of years and store as an attribute for later use in constructing
  ## projection matrices etc.
  attr(catch_df, "T") <- length(unique(catch_df$time))
  catch_df
}

##' Randomly choose a subset of observations for some or all vessels etc.
##'
##' Takes subsets of observations by year and other characteristics. Subsetting
##' is specified as a data frame. Any unspecified levels are not subset.
##' @title Subsample catch by year
##' @param catch_df Catch data frame to subset
##' @param n_df Data frame specifiying how to subset
##' @return A data frame with the same columns, subset as specified
##' @author John Best
##'
##' @examples
##' \dontrun{
##'   ## To subset the observations of a single vessel, only specify a number of
##'   ## observations for that vessel
##'   subsample_catch(catch_df,
##'                   tibble(vessel_idx = 2,
##'                          n = 1000))
##'   ## To subset the observations for two vessels, specify a number of
##'   ## observations for each
##'   subsample_catch(catch_df,
##'                   tibble(vessel_idx = c(1, 2),
##'                          n = c(200, 1000)))
##' }
##' @importFrom dplyr %>%
##' @export
subsample_catch <- function(catch_df, n_df = NULL) {
  if (is.null(n_df)) return(catch_df)
  nms <- names(n_df)
  "n" %in% names(n_df) || stop("Need column \`n\`")
  nms_join <- nms[nms != "n"]
  catch_df <- catch_df %>%
    tidyr::nest(data = c(-vessel_idx, -time)) %>%
    dplyr::left_join(n_df, by = nms_join) %>%
    dplyr::mutate(n = dplyr::coalesce(n, purrr::map_dbl(data, nrow))) %>%
    dplyr::mutate(data = purrr::map2(data, n, dplyr::sample_n)) %>%
    dplyr::select(-n) %>%
    tidyr::unnest(data)
  ## Double check that T is correct here
  attr(catch_df, "T") <- length(unique(catch_df$time))
  catch_df
}

##' Check that scenario name is correct. Available options are "pref", "spat",
##' or "combo". Any other string will throw an error.
##'
##' @title Check scenario name
##' @param sc Scenario name, as a string
##' @return \code{TRUE} if scenario is correct
##' @author John Best
##' @export
check_scenario <- function(sc) {
  sc %in% c("pref", "spat", "combo") ||
    stop("Scenario must be one of \"pref\", \"spat\", or \"combo\"")
  TRUE
}

##' Read true population total from file
##'
##' @title Read true population state for each year
##' @param repl Replicate number
##' @param sc Scenario; one of "pref", "spat", or "combo"
##' @param root_dir Directory with e.g. \code{repl_01} subdirectory that
##'   contains \code{popstate_*.csv}
##' @return A \code{tibble} with population and year, starting from 1
##' @author John Best
##' @export
read_popstate <- function(repl, sc, root_dir = ".") {
  check_scenario(sc)
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <- stringr::str_pad(repl, 2, pad = "0")
  repl_dir <- paste0("repl_", repl_str)
  sc_file <- paste0("popstate_", repl_str, "_", sc, ".csv")
  flnm <- file.path(root_dir, repl_dir, sc_file)
  pop <- readr::read_csv(flnm,
                         col_types = readr::cols(pop = readr::col_double()))
  dplyr::mutate(pop, time = seq_along(pop))
}

##' Assumes the simulation \eqn{[0, 100]×[0, 100]} domain.
##'
##' @title Define the domain boundary
##' @return An \code{inla.mesh.segment} giving the boundary of the 100 by 100
##'   domain
##' @author John Best
##' @export
domain_boundary <- function() {
  boundary <- INLA::inla.mesh.segment(matrix(c(0, 0,   100, 100, 0,
                                               0, 100, 100, 0,   0),
                                             ncol = 2))
}

##' Generate the coordinates of a regular grid on a \eqn{[0, 100]×[0, 100]}
##' domain with steps in both direction of \code{step}.
##'
##' @title Generate a location grid
##' @param step Distance between subsequent locations; defaults to 1
##' @return A two-column matrix with of locations in two dimensions, each
##'   \code{step} apart.
##' @author John Best
##' @export
loc_grid <- function(step = 1.0) {
  grid_start <- step / 2
  grid_end <- 100 - step / 2
  grid_seq <- seq(grid_start, grid_end, step)
  as.matrix(expand.grid(s1 = grid_seq,
                        s2 = grid_seq))
}

##' Generate the standard mesh used for simulation fits
##'
##' @title Generate INLA mesh
##' @return A \code{inla.mesh} object that can be passed to \code{generate_fem}.
##' @author John Best
##' @export
generate_mesh <- function() {
  ## Discretize spatial domain into mesh
  boundary <- domain_boundary()
  loc <- loc_grid(2.0)
  INLA::inla.mesh.2d(loc,
                     boundary = boundary,
                     offset = c(0.0, 30.0),
                     # Shortest correlation range should be ~30
                     max.edge = c(5, 20),
                     max.n = c(400, 100),
                     min.angle = c(30, 21),
                     cutoff = 5)
}

##' Generate the FEM matrices with appropriate names to pass to TMB as a
##' \code{spde_t}.
##'
##' @title Generate FEM matrices
##' @param mesh INLA mesh to generate FEM matrices from
##' @return List of sparse matrices suitable for passing to TMB.
##' @author John Best
##' @export
generate_fem <- function(mesh) {
  ## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended
  ## by Finn Lindgren
  fem <- INLA::inla.mesh.fem(mesh)
  ## Rename elements of `fem` so that they are recognized as a `spde` object in
  ## TMB
  names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")
  fem
}

##' Generate a projection matrix to locations in \code{data_df}, possible
##' zeroing out some observations. Also group by year for a spatiotemporal
##' effect.
##'
##' @title Generate a projection matrix
##' @param mesh INLA mesh to project
##' @param data_df Data frame with \code{s1} and \code{s2} columns as
##'   coordinates to project to
##' @param vessel_idx Vessel index to project to; others receive zero weight and
##'   are dropped
##' @param group Year of observation for spatiotemporal effect
##' @param zero Is this an empty projection matrix? (May be useful if e.g. the
##'   random effects parameters associated with it are map'd to zeros.)
##' @return A sparse projection matrix
##' @author John Best
##' @export
generate_projection <- function(mesh, data_df, vessel_idx = NULL, group = NULL,
                                zero = FALSE) {
  ## If effect is map'd, return an all-zero projection matrix. Should save some
  ## unnecessary multiplications?
  if (zero) return(generate_empty_projection(mesh, data_df, group))
  if (is.null(vessel_idx)) {
    wts <- NULL
  } else {
    wts <- data_df$vessel_idx %in% vessel_idx
  }
  loc <- as.matrix(data_df[, c("s1", "s2")])
  A <- INLA::inla.spde.make.A(mesh, loc, group = group, weights = wts)
  Matrix::drop0(A)
}

##' When spatial or spatiotemporal effects are map'd to zero, pass an empty
##' projection matrix. Not exported, use `zero = TRUE` in `generate_projection`
##' instead.
##'
##' @title Generate an empty project matrix
##' @param mesh INLA mesh to project
##' @param data_df Data frame with \code{s1} and \code{s2} columns as
##'   coordinates to project to
##' @param group Year of observation for spatiotemporal effect
##' @return Sparse \code{Matrix::Matrx} of the appropriate dimensions, but
##'   filled with zeros
##' @author John Best
generate_empty_projection <- function(mesh, data_df, group = NULL) {
  n_obs <- nrow(data_df)
  if (!is.null(group)) {
    n_years <- length(unique(group))
  } else {
    n_years <- 1
  }
  Matrix::Matrix(0, nrow = n_obs, ncol = mesh$n * n_years)
}

##' Parse coordinates that were saved as a tuple of \code{Float64} and read in
##' as a string. Much slower than using \code{get_coordref} and doing a
##' \code{join}.
##'
##' @title Parse a tuple of numbers read as a string
##' @param coord_tuple A string of the form \code{"(1.0, 2.0)"}
##' @return The two numbers in the tuple
##' @author John Best
##' @export
parse_coords <- function(coord_tuple) {
  coord_str <- stringr::str_sub(coord_tuple, start = 2L, end = -2L)
  coords <- as.numeric(stringr::str_split_fixed(coord_str, ", ", 2))
  coords
}

##' Get the table of location indices and corresponding location coordinates
##'
##' Assumes that the table is in the file \code{coordref.csv} in the working
##' directory, as saved by \code{FisherySim.jl}
##' @title Read the index to coordinate reference table
##' @param root_dir Directory containing the file "coordref.csv"
##' @return A \code{tibble} with columns \code{loc_idx}, \code{s1}, and
##'   \code{s2}.
##' @author John Best
##' @export
get_coordref <- function(root_dir = ".") {
  cr_file <- normalizePath(file.path(root_dir, "coordref.csv"), mustWork = TRUE)
  readr::read_csv(cr_file,
                  col_types = readr::cols(
                    loc_idx = readr::col_integer(),
                    s1 = readr::col_double(),
                    s2 = readr::col_double()))
}

##' Create a data frame with a regular grid over the spatial domain for each
##' year of the index.
##'
##' @title Create a data frame with index locations and times
##' @param step Step in the spatial domain
##' @param T Number of years
##' @return A \code{tibble} with columns \code{s1}, \code{s2}, and \code{time}
##'   and attribute \code{step} for calculating the integration weights.
##' @author John Best
##' @export
create_index_df <- function(step = 1.0, T = 25) {
  loc_df <- tibble::as_tibble(loc_grid(step)) %>%
   dplyr::mutate(loc_idx = seq_along(s1))
  index_df <- purrr::cross_df(list(loc_idx = loc_df$loc_idx,
                                   time = seq_len(T))) %>%
    dplyr::full_join(loc_df, by = "loc_idx") %>%
    dplyr::select(-loc_idx)
  attr(index_df, "step") <- step
  index_df
}

##' Create a list of all data elements required to fit a \code{spatq} model with
##' names etc.
##'
##' @title Prepare the \code{data} argument for the \code{spatq} model
##' @param catch_df A data frame with catch observations, as read in by
##'   \code{read_catch}
##' @param index_df A data frame with coordinates for index calculation, as
##'   generated by \code{create_index_df}
##' @param mesh A FEM mesh for the domain
##' @param fem The sparse FEM matrices from the \code{mesh} as generated by
##'   \code{generate_fem}
##' @return A \code{list} with the data values for a \code{spatq} model
##' @author John Best
##' @export
prepare_data <- function(catch_df, index_df, mesh, fem) {
  T <- length(unique(catch_df$time))

  ## If single-vessel and no covariates, nothing to estimate for `lambda_n` or
  ## `lambda_w`. Check length of unique column contents rather than `levels` in
  ## case already a factor and has all original vessels as levels even if none
  ## present. Useful e.g. if estimating index from single survey vessel.
  if (length(unique(catch_df$vessel_idx)) >= 2) {
    ## Catchability fixed effect design matrix. Drop first column so
    ## that fixed effects are identifiable (e.g. if only two vessels,
    ## only adjust catchability for second.) Use `drop = FALSE` so the
    ## result is a matrix if only one column is left.
    R_n <- stats::model.matrix(~ factor(vessel_idx),
                               data = catch_df)[, -1, drop = FALSE]
    R_w <- stats::model.matrix(~ factor(vessel_idx),
                               data = catch_df)[, -1, drop = FALSE]
  } else {
    R_n <- matrix(0, nrow = nrow(catch_df))
    R_w <- matrix(0, nrow = nrow(catch_df))
  }
  dat <- list(catch_obs = catch_df$catch_biomass,
              area_swept = catch_df$effort,
              ## Abundance fixed effect design matrices
              X_n = stats::model.matrix(~ factor(time) + 0, data = catch_df),
              X_w = stats::model.matrix(~ factor(time) + 0, data = catch_df),
              IX_n = stats::model.matrix(~ factor(time) + 0, data = index_df),
              IX_w = stats::model.matrix(~ factor(time) + 0, data = index_df),
              ## Abundance random effect design matrices
              Z_n = matrix(0, nrow = nrow(catch_df), ncol = 1),
              Z_w = matrix(0, nrow = nrow(catch_df), ncol = 1),
              IZ_n = matrix(0, nrow = nrow(index_df), ncol = 1),
              IZ_w = matrix(0, nrow = nrow(index_df), ncol = 1),
              ## Abundance projection matrices
              A_spat = generate_projection(mesh, catch_df),
              A_sptemp = generate_projection(mesh, catch_df,
                                             group = catch_df$time),
              IA_spat = generate_projection(mesh, index_df),
              IA_sptemp = generate_projection(mesh, index_df,
                                              group = index_df$time),
              ## Integration weights
              Ih = rep_len(attr(index_df, "step")^2, nrow(index_df)),

              ## Catchability fixed effect design matrix. See above for
              ## specification.
              R_n = R_n,
              R_w = R_w,
              ## Catchability random effect design matrix
              V_n = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              V_w = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              ## Catchability projection matrices
              A_qspat = generate_projection(mesh, catch_df, vessel_idx = 2),
              A_qsptemp = generate_projection(mesh, catch_df,
                                              vessel_idx = 2,
                                              group = catch_df$time,
                                              zero = TRUE),

              ## FEM matrices for spatial/spatiotemporal effects
              spde = fem,

              ## Turn off random effects likelihood components. Start with all
              ## on, turn off after map
              proc_switch = rep(TRUE, 6),
              ## Normalize GMRFs externally, and include flag for returning
              ## negative log-likelihood without data likelihood
              norm_flag = FALSE,
              incl_data = TRUE)
  ## Store the number of years as an attribute
  attr(dat, "T") <- T
  dat
}

##' Generate an (all-zero) parameter starting value of appropriate dimension
##' from the corresponding data element
##'
##' @title Parameter values from data elements
##' @param data_el Data element to take dimensions from
##' @param T Number of years
##' @return A number, vector, or matrix of zeros of appropriate dimension and
##'   type for that parameter of the \code{spatq} model
##' @author John Best
##' @export
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

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. The kappa parameter is \code{sqrt(8) / rho}.
##'
##' @title Calculate kappa from rho
##' @param rho Correlation range parameter
##' @return The value of kappa corresponding to rho
##' @author John Best
##' @export
pars_kappa <- function(rho) {
  sqrt(8) / rho
}

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. The tau parameter is \code{sqrt(1 / (sig2 * 4
##' * pi * kappa^2))}.
##'
##' @title Calculate tau from rho and sigma^2
##' @param sig2 Marginal variance
##' @param rho Correlation range parameter
##' @return The value of tau corresponding to the specified correlation range
##'   and marginal variance
##' @author John Best
##' @export
pars_tau <- function(sig2, rho) {
  den <- sig2 * 4 * pi * pars_kappa(rho)^2
  sqrt(1 / den)
}

##' Estimates numbers density fixed effects through a binomial regression with
##' the complementary log-log link function. Log numbers density and log
##' encounter probability provide an offset for a log normal regression on
##' positive catches, hopefully providing reasonable starting parameter values
##' for weight per group covariates. Log numbers density is NOT currently
##' bias-corrected, but should still work well as starting values.
##'
##' @title Initial estimates of fixed effect parameters
##' @param data A data list, as produced by \code{prepare_data}
##' @return A list with estimates of \code{beta_n}, \code{beta_w},
##'   \code{lambda_n}, and \code{lambda_w}
##' @author John Best
##' @export
init_fixef <- function(data) {
  ## Don't include catchability effects if they aren't used in the model (e.g.
  ## if estimating index using single survey vessel).
  q_used <- any(data$R_n != 0)
  if (q_used) {
    df_n <- tibble::tibble(enc = data$catch_obs > 0,
                           X_n = data$X_n,
                           R_n = data$R_n)
  } else {
    df_n <- tibble::tibble(enc = data$catch_obs > 0,
                           X_n = data$X_n)
  }
  ## Estimate each group density fixed effect with a simple GLM, using the
  ## complementary log-log link for encounter. This corresponds directly to
  ## group density in the Poisson link model.
  mod_n <- stats::glm(enc ~ 0 + ., data = df_n,
                      family = stats::binomial(cloglog))
  est_n <- stats::coef(mod_n)

  ## Using the estimated group density and probability of encounter, calculate
  ## the log mean weight per group, conditional on a positive catch.
  enc <- df_n$enc
  if (q_used) {
    df_w <- tibble::tibble(catch_obs = data$catch_obs[enc],
                           X_w = data$X_w[enc, ],
                           R_w = data$R_w[enc, ])
  } else {
    df_w <- tibble::tibble(catch_obs = data$catch_obs[enc],
                           X_w = data$X_w[enc, ])
  }
  ## Log-positive catch rate is log(n) - log(p) + log(w), so the first two are
  ## used as an offset here.
  offset_w <- stats::predict(mod_n) -
    log(stats::predict(mod_n, type = "response"))
  offset_w <- offset_w[data$catch_obs > 0]
  ## mod_w <- glm(catch_obs ~ 0 + ., data = df_w, offset = offset_w,
  ##              family = gaussian(log))
  mod_w <- stats::lm(log(catch_obs) ~ 0 + ., data = df_w, offset = offset_w)
  ## Not going to worry about the bias correction for these initial values
  est_w <- stats::coef(mod_w)

  init <- list(beta_n = est_n[seq_len(ncol(data$X_n))],
               beta_w = est_w[seq_len(ncol(data$X_w))],
               lambda_n = ifelse(q_used, est_n[-seq_len(ncol(data$X_n))], 0),
               lambda_w = ifelse(q_used, est_w[-seq_len(ncol(data$X_w))], 0))
  lapply(init, unname)
}

##' Prepare a list of parameters of appropriate dimension. Otherwise
##' uninitialized parameter starting values are set to zeros except kappa and
##' tau, which use a correlation range of 50 (half the domain) and marginal
##' variance of 1, tranformed to the kappa and tau parameterization. Initial
##' estimates of fixed effects can be estimated using GLMs with the
##' \code{init_fixef = TRUE} argument. This is recommended.
##'
##' @title Prepare parameter list
##' @param data Data list, as from \code{prepare_data}
##' @param mesh SPDE mesh
##' @param init_fixef Use GLMs to estimate fixed effect parameters for abundance
##'   and catchability in both numbers density and weight per group processes
##' @return List with scalar, vector, or matrices of zeros as starting values
##'   for optimization
##' @author John Best
##' @export
prepare_pars <- function(data, mesh, init_fixef = TRUE) {
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
  if (init_fixef) {
    init_est <- init_fixef(data)
    pars$beta_n <- init_est$beta_n
    pars$beta_w <- init_est$beta_w
    pars$lambda_n <- init_est$lambda_n
    pars$lambda_w <- init_est$lambda_w
  }
  ## Add attribute noting whether `lambda_n` and `lambda_w` are map'd or not.
  ## Useful for consistency checks later.
  attr(pars, "map_lambda") <- all(data$R_n == 0)
  pars
}

##' Spatial/spatiotemporal parameters (kappa and tau) are stored in a length-8
##' vector. This function takes a name and returns the corresponding index.
##'
##' @title Spatial parameter index
##' @param par_name Name of the spatial/spatiotemporal parameter (e.g.
##'   \code{omega_n}, \code{phi_w} &c.)
##' @return The index of the effect in the parameter vectors
##' @author John Best
##' @export
spat_par_idx <- function(par_name) {
  spat_pars <- c("omega_n", "omega_w",
                 "epsilon_n", "epsilon_w",
                 "phi_n", "phi_w",
                 "psi_n", "psi_w")
  which(par_name == spat_pars)
}

##' The variance of general random effects are stored in a length-4 vector. This
##' function takes a name and returns the corresponding index.
##'
##' @title General random effects parameter index
##' @param par_name Name of the spatial/spatiotemporal parameter (e.g.
##'   \code{gamm_n}, \code{eta_w} &c.)
##' @return The index of the effect in the parameter vectors
##' @author John Best
##' @export
re_par_idx <- function(par_name) {
  re_pars <- c("gamma_n", "gamma_w",
               "eta_n", "eta_w")
  which(par_name == re_pars)
}

##' Prepare a \code{map} argument for \code{MakeADFun}.
##'
##' @title Prepare a \code{map}
##' @param pars Parameter list, as from \code{prepare_pars}
##' @param spec List of logicals indicating which parameters are to be
##'   estimated, as output by \code{\link{specify_estimated}}
##' @return A \code{map} list suitable for passing to \code{MakeADFun}
##' @author John Best
##' @export
prepare_map <- function(pars, spec) {
  ## Drop parameters that are *not* mapped
  spec2 <- spec[!vapply(spec, all, TRUE)]
  ## Invert to specify *mapped* parameter vectors of correct length
  mapped <- lapply(names(spec2),
                   function(nm) {
                     rep(!spec2[[nm]], length.out = length(pars[[nm]]))
                   })
  names(mapped) <- names(spec2)
  map <- lapply(mapped,
                function(mpd) {
                  ## Add individual levels for unmapped parameters
                  v <- cumsum(!mpd)
                  ## Specify NA for mapped parameters
                  v[mpd] <- NA
                  ## Convert to factor vector for TMB
                  factor(v)
                })
  return(map)
}

##' Prepare a character vector indicating which parameters should be
##' marginalized out using the Laplace approximation. Does not include
##' parameters that are \code{map}'d.
##'
##' Parameters in this model that are typically marginalized are:
##' \itemize{
##'   \item \code{gamma_n}, \code{gamma_w}
##'   \item \code{omega_n}, \code{omega_w}
##'   \item \code{epsilon_n}, \code{epsilon_w}
##'   \item \code{eta_n}, \code{eta_w}
##'   \item \code{phi_n}, \code{phi_w}
##'   \item \code{psi_n}, \code{psi_w}
##' }
##'
##' @title Prepare \code{random}
##' @param map A \code{map} list, as from \code{prepare_map}
##' @return A character vector indicating which parameters should be integrated
##'   out using the Laplace approximation.
##' @author John Best
##' @export
prepare_random <- function(map) {
  re_pars <- c("gamma_n", "gamma_w",
               "omega_n", "omega_w",
               "epsilon_n", "epsilon_w",
               "eta_n", "eta_w",
               "phi_n", "phi_w",
               "psi_n", "psi_w")
  ## Check that all `map` random effects parameter entries are all NAs; check
  ## that none are included but not actually map'd
  if (!all(vapply(re_pars, function(p) all(is.na(map[[p]])),
                  FUN.VALUE = TRUE)))
    stop("Map'd random effects parameters must be all NAs")
  ## Return vector of random effects names with map'd parameter names removed
  setdiff(re_pars, names(map))
}

##' Generates a vector of length 6 indicating whether each pair of numbers
##' density and weight per group random processes should be included in the
##' likelihood calculation; i.e. they are map'd.
##'
##' @title Prepare process switch
##' @param random Character vector of random parameters; as from
##'   \code{prepare_random}
##' @return Logical vector of length 6 indicating which random processes are
##'   *not* map'd off
##' @author John Best
prepare_proc_switch <- function(random) {
  procs <- c("gamma_n", "gamma_w", "omega_n", "omega_w",
             "epsilon_n", "epsilon_w", "eta_n", "eta_w",
             "phi_n", "phi_w", "psi_n", "psi_w")
  on <- procs %in% random
  vapply(seq(1, 11, 2), function(i) on[i] || on[i + 1], TRUE)
}

##' Verify data, parameters, and map, then contruct the ADFun.
##'
##' @title Make a TMB ADFun
##' @param data Data list, as from \code{prepare_data}
##' @param parameters Parameter list, as from \code{prepare_pars}
##' @param map A \code{map} list, as from \code{prepare_map}
##' @param random A character vector of parameters to be marginalized, as from
##'   \code{prepare_random}
##' @param ... Additional arguments to pass to \code{MakeADFun}
##' @param DLL TMB model to use; only used for debugging
##' @param silent Output TMB progress?
##' @param runSymbolicAnalysis Use Metis reorderings? (Requires special
##'   installation of TMB; see documentation.)
##' @param normalize Normalize GMRF likelihoods?
##' @return A TMB ADFun suitable for optimization
##' @author John Best
##' @export
prepare_adfun <- function(data, parameters, map, random,
                          ..., DLL = "spatq", silent = TRUE,
                          runSymbolicAnalysis = TRUE,
                          normalize = FALSE) {
  data$proc_switch <- prepare_proc_switch(random)
  data$norm_flag <- normalize
  verify_spatq(data, parameters, map)
  obj <- TMB::MakeADFun(data = data,
                        parameters = parameters,
                        map = map,
                        random = random,
                        DLL = DLL,
                        ...,
                        silent = silent)
  if (!normalize)
    obj <- TMB::normalize(obj, flag = "incl_data", value = FALSE)
  if (runSymbolicAnalysis)
    TMB::runSymbolicAnalysis(obj)
  obj
}

##' Read in a simulated data set and construct a TMB ADFun for model fitting.
##' From the working directory, replicates should be in subdirectories and
##' scenarios named appropriately. See \code{read_catch} for details.
##'
##' @title Create ADFun from simulated data
##' @param repl Replicate number
##' @param sc Scenario; "pref", "spat", or "combo"
##' @param sub_df Data frame indicating subsampling strategy; see
##'   \code{subsample_catch}
##' @param root_dir Directory to load data from
##' @param max_T Last year of data to include
##' @param spec_estd List of logicals indicating which parameters are to be
##'   estimated, as output \code{\link{specify_estimated}}
##' @param ... Additional options to pass to \code{prepare_adfun}
##' @return A TMB ADFun suitable for optimization
##' @author John Best
##' @export
make_sim_adfun <- function(repl, sc, sub_df = NULL,
                           root_dir = ".", max_T = NULL,
                           spec_estd = specify_estimated(), ...) {
  ## Read in data
  catch_df <- read_catch(repl, sc, root_dir)
  if (!is.null(max_T)) {
    catch_df <- dplyr::filter(catch_df, time <= max_T)
    attr(catch_df, "T") <- max_T
  }
  ## Subset observations
  catch_df <- subsample_catch(catch_df, sub_df)

  ## Create index integration reference
  index_df <- create_index_df(step = 5, T = attr(catch_df, "T"))

  ## Discretize space
  mesh <- generate_mesh()
  fem <- generate_fem(mesh)

  ## Prepare model specification components
  data <- prepare_data(catch_df, index_df, mesh, fem)
  parameters <- prepare_pars(data, mesh)
  map <- prepare_map(parameters,
                     spec = spec_estd)
  random <- prepare_random(map)

  ## Verify and construct ADFun
  prepare_adfun(data, parameters, map, random, ...)
}
