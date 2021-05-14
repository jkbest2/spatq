##' Read a CSV containing the true and estimated indices of abundance.
##'
##' Includes true, estimated, and debiased estimates as well as error estimates.
##' @title Read indices of abundance
##' @param filestudy Either the file path or the study name
##' @param repl Replicate number
##' @param opmod Operating model number
##' @param estmod Estimation model name
##' @param root_dir Root directory
##' @param estmods Names of estimation models to code column as a factor
##' @return \code{\link[tibble]{tibble}} with true and estimated indices of
##'   abundance
##' @author John K Best
##' @export
read_index_csv <- function(filestudy,
                           repl = NULL,
                           opmod = NULL,
                           estmod = NULL,
                           root_dir = ".",
                           estmods = c("survey",
                                       "spatial_ab",
                                       "spatial_q")) {
  ## Construct filename if only components are specified
  if (!is.null(repl) && !is.null(opmod) && !is.null(estmod)) {
    filename <- res_file_paths(filestudy,
                               repl,
                               opmod,
                               estmod,
                               root_dir)$indexcsv
  } else {
    filename <- filestudy
  }

  readr::read_csv(filename,
                  col_types = readr::cols(
                    study = readr::col_character(),
                    repl = readr::col_integer(),
                    opmod = readr::col_integer(),
                    estmod = readr::col_factor(levels = estmods),
                    year = readr::col_integer(),
                    raw_est = readr::col_double(),
                    index_est = readr::col_double(),
                    raw_unb = readr::col_double(),
                    index_unb = readr::col_double(),
                    raw_sd = readr::col_double(),
                    index_sd = readr::col_double(),
                    raw_unb_sd = readr::col_double(),
                    unb_sd = readr::col_double(),
                    raw_true = readr::col_double(),
                    index_true = readr::col_double())) %>%
    dplyr::mutate(repl = factor(repl),
                  opmod = factor(opmod))
}
##' Read all fitted index CSVs
##'
##' @title Load index CSVs for a study
##' @param study Simulation study name
##' @param repls Replicate range
##' @param opmods Operating model range
##' @param estmods Estimation model names
##' @param root_dir Root directory
##' @return \code{\link[tibble]{tibble}} with indices for all fitted models
##' @author John K Best
##' @export
read_all_indices <- function(study,
                             repls,
                             opmods,
                             estmods,
                             root_dir = ".") {
  purrr::cross(list(study = study,
                    repl = repls,
                    opmod = opmods,
                    estmod = estmods,
                    root_dir = root_dir)) %>%
    purrr::map_df(read_index, filetype = "csv")
}

##' @title Read an index from a fitted model
##' @param x Path to index or \code{\link{spatq_simstudyspec}}
##' @param filetype File type, either "feather" or "csv"
##' @param estmods Names of estimation models
##' @return A data frame with fitted index values
##' @author John K Best
##' @export
read_index <- function(x, filetype, estmods) UseMethod("read_index")
##' @export
read_index.character <- function(x,
                                 filetype = NULL,
                                 estmods = c("survey", "spatial_ab", "spatial_q")) {
  if (is.null(filetype)) {
    if (grepl("\\.feather$", x)) {
      filetype <- "feather"
    } else if (grepl("\\.csv$")) {
      filetype <- "csv"
    }
  }

  if (filetype == "feather") {
    index_df <- arrow::read_feather(x) %>%
      dplyr::mutate(estmod = factor(estmod, levels = estmods))
  } else if (filetype == "csv") {
    index_df <- readr::read_csv(x,
                                col_types =
                                  readr::cols(
                                           study = readr::col_character(),
                                           repl = readr::col_integer(),
                                           opmod = readr::col_integer(),
                                           estmod = readr::col_factor(levels = estmods),
                                           year = readr::col_integer(),
                                           raw_est = readr::col_double(),
                                           index_est = readr::col_double(),
                                           raw_unb = readr::col_double(),
                                           index_unb = readr::col_double(),
                                           raw_sd = readr::col_double(),
                                           index_sd = readr::col_double(),
                                           raw_unb_sd = readr::col_double(),
                                           unb_sd = readr::col_double(),
                                           raw_true = readr::col_double(),
                                           index_true = readr::col_double()))
  }
  index_df %>%
    dplyr::mutate(repl = factor(repl),
                  opmod = factor(opmod))
}
##' @export
read_index.spatq_simstudyspec <- function(x,
                                          filetype = "csv",
                                          estmods = c("survey", "spatial_ab", "spatial_q")) {
  file <- index_path(x, filetype)
  read_index(file, filetype, estmods)
}
##' @export
read_index.list <- function(x,
                            filetype = "csv",
                            estmods = c("survey", "spatial_ab", "spatial_q")) {
  spec <- spatq_simstudyspec(x)
  read_index(spec, filetype, estmods)
}

##' @title Read saved Rdata from a fitted model
##' @param x Path to index, \code{\link{spatq_simstudyspec}}, or list
##' @return A data frame with fitted index values
##' @author John K Best
##' @export
read_rdata <- function(x) UseMethod("read_rdata")
##' @export
read_rdata.character <- function(x) {
  readRDS(x)
}
##' @export
read_rdata.spatq_simstudyspec <- function(x) {
  path <- rdata_path(x)
  read_rdata(path)
}
##' @export
read_rdata.list <- function(x) {
  path <- rdata_path(x)
  read_rdata(path)
}
