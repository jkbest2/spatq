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
read_all_indices <- function (study,
                              repls,
                              opmods,
                              estmods,
                              root_dir = ".") {
  csv_list <- all_res_file_paths(study, repls, opmods, estmods, root_dir)$indexcsv
  purrr::map_df(csv_list, read_index_csv)
}
