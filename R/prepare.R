##' Read in a replicate of data simulated by FisherySim.jl
##'
##' Requires that catch files are e.g. \code{repl_01/catch_simple.csv}.
##' @title Read simulated data
##' @param repl Replicate number
##' @param sc Scenario; one of "naive", "simple", "scaled", or "shared"
##' @param dir Directory to work in; defaults to current working directory
##' @return A \code{tibble} with catch observations
##' @author John Best
##' @export
read_catch <- function(repl, sc, dir = ".") {
  sc %in% c("naive", "simple", "scaled", "shared")
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <- str_pad(repl, 2, pad = "0")
  flnm <- paste0(dir, "repl_", repl_str, "/catch_", repl_str, "_", sc, ".csv")
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
##' @export
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

##' Read true population total from file
##'
##' @title Read true population state for each year
##' @param repl Replicate number
##' @param sc Scenario; one of "naive", "simple", "scaled", or "shared"
##' @param dir Directory to start from; defaults to current directory
##' @return A \code{tibble} with population and year, starting from 1
##' @author John Best
##' @export
read_popstate <- function(repl, sc) {
  sc %in% c("naive", "simple", "scaled", "shared")
  ## File names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two
  ## digits.
  repl_str <-str_pad(repl, 2, pad = "0")
  flnm <- paste0("repl_", repl_str, "/popstate_", repl_str, "_", sc, ".csv")
  pop <-read_csv(flnm, col_types =cols(pop =col_double()))
 mutate(pop, time = seq_along(pop))
}

##' Assumes the simulation [0, 100] x [0, 100] domain.
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

##' Generate the coordinates of a regular grid on a [0, 100] x [0, 100] domain
##' with steps in both direction of \code{step}.
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
                     offset = c(5.0, 25.0),
                     max.edge = c(3 * sqrt(2), 6 * sqrt(2)),
                     min.angle = c(30, 21),
                     cutoff = 0.5)
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
##' @return A sparse projection matrix
##' @author John Best
##' @export
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

##' When spatial or spatiotemporal effects are map'd to zero, pass an empty
##' projection matrix.
##'
##' @title Generate an empty project matrix
##' @param mesh INLA mesh to project
##' @param data_df Data frame with \code{s1} and \code{s2} columns as
##'   coordinates to project to
##' @param group Year of observation for spatiotemporal effect
##' @param ... Not used but availble for consistent calls with
##'   \code{generate_projection}
##' @return Sparse \code{Matrix::Matrx} of the appropriate dimensions, but
##'   filled with zeros
##' @author John Best
##' @export
generate_empty_projection <- function(mesh, data_df, group = NULL, ...) {
  n_obs <- nrow(data_df)
  if (!is.null(group)) {
    n_years <- length(unique(group))
  } else {
    n_years = 1
  }
  Matrix::Matrix(0, nrow = mesh$n, ncol = n_obs * n_years)
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
  coord_str <- str_sub(coord_tuple, start = 2L, end = -2L)
  coords <- as.numeric(str_split_fixed(coord_str, ", ", 2))
  coords
}

##' Get the table of location indices and corresponding location coordinates
##'
##' Assumes that the table is in the file \code{coordref.csv} in the working
##' directory, as saved by \code{FisherySim.jl}
##' @title Read the index to coordinate reference table
##' @return A \code{tibble} with columns \code{loc_idx}, \code{s1}, and
##'   \code{s2}.
##' @author John Best
##' @export
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
  loc_df <- as_tibble(loc_grid(step)) %>%
   mutate(loc_idx = seq_along(s1))
  index_df <- cross_df(list(loc_idx = loc_df$loc_idx, time = seq_len(T))) %>%
    full_join(loc_df, by = "loc_idx") %>%
    select(-loc_idx)
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
  dat <- list(catch_obs = catch_df$catch_biomass,
              area_swept = catch_df$effort,
              ## Abundance fixed effect design matrices
              X_n = model.matrix(~ factor(time) + 0, data = catch_df),
              X_w = model.matrix(~ factor(time) + 0, data = catch_df),
              IX_n = model.matrix(~ factor(time) + 0, data = index_df),
              IX_w = model.matrix(~ factor(time) + 0, data = index_df),
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

              ## Catchability fixed effect design matrix
              R_n = model.matrix(~ factor(vessel_idx), data = catch_df),
              R_w = model.matrix(~ factor(vessel_idx), data = catch_df),
              ## Catchability random effect design matrix
              V_n = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              V_w = matrix(0.0, nrow = nrow(catch_df), ncol = 1),
              ## Catchability projection matrices
              A_qspat = generate_projection(mesh, catch_df, vessel_idx = 2),
              A_qsptemp = generate_projection(mesh, catch_df,
                                            vessel_idx = 2,
                                            group = catch_df$time),
              ## FEM matrices for spatial/spatiotemporal effects
              spde = fem)
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

##' Prepare a list of parameters of appropriate dimension. All starting values
##' are zeros except kappa and tau, which use a correlation range of 50 (half
##' the domain) and marginal variance of 1, tranformed to the kappa and tau
##' parameterization.
##'
##' @title Prepare parameter list
##' @param data Data list, as from \code{prepare_data}
##' @param mesh SPDE mesh
##' @return List with scalar, vector, or matrices of zeros as starting values
##'   for optimization
##' @author John Best
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

##' Prepare a \code{map} argument for \code{MakeADFun}. Accepts parameter names
##' and takes care of dependencies such as spatial or random effect
##' hyperparameters being map'd when their effects are.
##'
##' @title Prepare a \code{map}
##' @param pars Parameter list, as from \code{prepare_pars}
##' @param map_pars Character vector indicating which parameters to \code{map}
##' @return A \code{map} list suitable for passing to \code{MakeADFun}
##' @author John Best
##' @export
prepare_map <- function(pars, map_pars) {
  map <- lapply(pars, function(par) par[] <- factor(seq_along(par)))
  for (par in map_pars) {
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

##' Prepare a character vector indicating which parameters should be
##' marginalized out using the Laplace approximation. Does not include
##' parameters that are \code{map}d.
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
  is_mapd <- unlist(lapply(re_pars, function(par) all(is.na(map[[par]]))))
  re_pars[!is_mapd]
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
##' @param silent Output TMB progress?
##' @return A TMB ADFun suitable for optimization
##' @author John Best
##' @export
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

##' Read in a simulated data set and construct a TMB ADFun for model fitting.
##' From the working directory, replicates should be in subdirectories and
##' scenarios named appropriately. See \code{read_catch} for details.
##'
##' @title Create ADFun from simulated data
##' @param repl Replicate number
##' @param sc Scenario; "naive", "simple", "scaled", or "shared"
##' @param sub_df Data frame indicating subsampling strategy; see
##'   \code{subsample_catch}
##' @param map A \code{map} list, as from \code{prepare_map}
##' @return A TMB ADFun suitable for optimization
##' @author John Best
##' @export
make_sim_adfun <- function(repl, sc, sub_df = NULL) {
  ## Read in data
  catch_df <- read_catch(repl, sc)
  ## Subset observations
  catch_df <- subsample_catch(catch_df, sub_df)

  ## Create index integration reference
  index_df <- create_index_df(step = 5, T = 25)

  ## Discretize space
  mesh <- generate_mesh()
  fem <- generate_fem(mesh)

  ## Prepare model specification components
  data <- prepare_data(catch_df, index_df, mesh, fem)
  parameters <- prepare_pars(data, mesh)
  map <- prepare_map(parameters,
                     map_list = c("gamma_n", "gamma_w",
                                  "epsilon_n", "epsilon_w",
                                  "eta_n", "eta_w",
                                  "psi_n", "psi_w"))
  random <- prepare_random(map)

  ## Verify and construct ADFun
  prepare_adfun(data, parameters, map, random)
}
