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
                                coord1 = readr::col_double(),
                                coord2 = readr::col_double(),
                                effort = readr::col_double(),
                                catch_biomass = readr::col_double())) %>%
    dplyr::rename(s1 = coord1, s2 = coord2)
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
##'   \code{\link{read_catch}}
##' @param index_df A data frame with coordinates for index calculation, as
##'   generated by \code{\link{create_index_df}}
##' @param mesh A FEM mesh for the domain
##' @param fem The sparse FEM matrices from the \code{mesh} as generated by
##'   \code{generate_fem}
##' @param obs_lik Observation likelihood. Current options are \code{0}:
##'   Poisson-link log-normal and \code{1}: Tweedie
##' @param X_contr Contrasts to use for fixed effect design matrix
##' @return A \code{list} with the data values for a \code{spatq} model
##' @author John Best
##' @importFrom stats contr.helmert
##' @export
prepare_data <- function(catch_df, index_df, mesh, fem,
                         obs_lik = 0L,
                         X_contr = contr.helmert) {
  T <- length(unique(catch_df$time))
  ## Convert `time` to a factor so different contrasts can be set easily. Don't
  ## overwrite `time` because later calls to `generate_projection` need it to be
  ## numeric.
  catch_df$ftime <- factor(catch_df$time)
  index_df$ftime <- factor(index_df$time)

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
              ## Abundance fixed effect design matrices, with `time` as a factor
              ## (`ftime`)
              X_n = stats::model.matrix(~ ftime, data = catch_df,
                                        contrasts.arg = list(ftime = X_contr)),
              X_w = stats::model.matrix(~ ftime, data = catch_df,
                                        contrasts.arg = list(ftime = X_contr)),
              IX_n = stats::model.matrix(~ ftime, data = index_df,
                                         contrasts.arg = list(ftime = X_contr)),
              IX_w = stats::model.matrix(~ ftime, data = index_df,
                                         contrasts.arg = list(ftime = X_contr)),
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
                                              group = catch_df$time),

              ## FEM matrices for spatial/spatiotemporal effects
              spde = fem,

              ## Turn off random effects likelihood components. Start with all
              ## on, turn off after map
              proc_switch = rep(TRUE, 12),

              ## Choose observation likelihood
              ## 0: Poisson-link zero-inflated log-normal
              ## 1: Tweedie
              obs_lik = obs_lik,

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
    par_el <- rep_len(0.0, n * T)
  } else if (is.null(dim(data_el))) {
    par_el <- 0.0
  } else {
    par_el <- rep_len(0.0, ncol(data_el))
  }
  par_el
}

##' Generates a vector of length 6 indicating whether each pair of numbers
##' density and weight per group random processes should be included in the
##' likelihood calculation; i.e. they are map'd.
##'
##' @title Prepare process switch
##' @param random Character vector of random parameters; as from
##'   \code{prepare_random}
##' @return Logical vector of length 12 indicating which random processes are
##'   *not* map'd off
##' @author John Best
prepare_proc_switch <- function(random) {
  procs <- c("gamma_n", "gamma_w", "omega_n", "omega_w",
             "epsilon_n", "epsilon_w", "eta_n", "eta_w",
             "phi_n", "phi_w", "psi_n", "psi_w")
  procs %in% random
}
