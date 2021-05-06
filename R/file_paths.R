##' Each set of simulations is stored in a separate directory; this function
##' serves to look up the name of that directory based on the name of that
##' directory.
##'
##' @title Directory where a study's simulations live
##' @param study Study name
##' @param root_dir Where to start
##' @return Directory path or file base name for study
##' @author John K Best
study_dir <- function(study, root_dir = ".") {
  sdir <- switch(study,
                 qdevscaling = "qdevscaling",
                 sharedq = "sharedq",
                 prefintensity = "prefintensity",
                 densdepq = "densdepq",
                 counterpref = "counterpref")
  file.path(root_dir, sdir)
}

##' @describeIn study_dir The file base for each study
study_file_base <- function(study) {
  study_file <- switch(study,
                       qdevscaling = "qdevscale_",
                       sharedq = "sharedq_",
                       prefintensity = "prefintensity_",
                       densdepq = "densdepq_",
                       counterpref = "counterpref_")
  study_file
}

##' Convenience function for writing replicate directory names
##'
##' @title Directory name for each replicate
##' @param repl Replicate number
##' @return A string with the number left-padded to two characters and
##'   concatenated with \code{"repl_"}
##' @author John K Best
repl_dir <- function(repl) {
  repl <- stringr::str_pad(repl, 2, pad = 0)
  paste0("repl_", repl)
}

##' @describeIn sim_file_paths The name of each simulation file
sim_file_names <- function(study, opmod) {
  study_file <- study_file_base(study)
  opmod <- stringr::str_pad(opmod, 2, pad = "0")
  list(catch_csv = paste0(study_file, opmod, "_catch.csv"),
       catch_feather = paste0(study_file, opmod, "_catch.feather"),
       pop_csv = paste0(study_file, opmod, "_popstate.csv"),
       pop_feather = paste0(study_file, opmod, "_popstate.feather"),
       pop_h5 = paste0(study_file, opmod, "_popstate.h5"))
}

##' Get the paths to each simulation file.
##
##' File names are ${study_dir}/repl_${repl}/${study_opmod}_catch.csv File
##' names are repl_$repl/catch_$repl_$sc.csv, with $repl padded to two digits.
##' @title Paths to simulated data sets
##' @param study Study
##' @param repl Replicate number
##' @param opmod Operating model
##' @param root_dir Where to start
##' @return A list of file paths with elements \code{catch}, \code{popcsv}, and
##'   \code{poph5}
##' @author John K Best
##' @export
sim_file_paths <- function(study, repl, opmod, root_dir = ".") {
  files <- sim_file_names(study, opmod)
  paths <- file.path(study_dir(study, root_dir),
                     repl_dir(repl),
                     files)
  names(paths) <- names(files)
  as.list(paths)
}

##' @describeIn res_file_paths Create a results directory if necessary
##' @export
create_res_dir <- function(study, repl = NULL, root_dir = ".") {
  res <- file.path(study_dir(study, root_dir), "results")
  if (!dir.exists(res))
    dir.create(res)
  if (!is.null(repl)) {
    res_repl <- file.path(study_dir(study, root_dir),
                           "results",
                           repl_dir(repl))
    ## If ${study_dir}/resutls/${repl_dir} doesn't exist, create it
    vapply(res_repl,
           function(rdir) {
             v <- TRUE
             if (!dir.exists(rdir)) {
               v <- dir.create(rdir)
             }
             v
           }, TRUE)
  }
  invisible(TRUE)
}

##' @describeIn res_file_paths File names for each study, operating, and
##'   estimation model combination
res_file_names <- function(study, opmod, estmod) {
  opmod <- stringr::str_pad(opmod, 2, pad = 0)
  res_files <- paste0(study_file_base(study),
                      opmod, "_",
                      estmod,
                      c(".Rdata", "_index.csv"))
  names(res_files) <- c("rdata", "indexcsv")
  res_files
}

##' Get paths for the Rdata and index CSV files resulting from a model fit. For
##' \code{all_res_file_paths}, paths will be returned sorted by replicate number
##' then operating model, then estimation model.
##'
##' @title Paths to each results file
##' @param study Study name
##' @param repl Replicate number(s)
##' @param opmod Operating model
##' @param estmod Estimation model
##' @param root_dir Root directory
##' @return A list with elements \code{rdata} and \code{indexcsv} with relative
##'   paths starting from \code{root_dir}
##' @author John K Best
##' @export
res_file_paths <- function(study, repl, opmod, estmod, root_dir = ".") {
  res_dir <- file.path(study_dir(study, root_dir),
                       "results",
                       repl_dir(repl))
  res_files <- res_file_names(study, opmod, estmod)
  res_paths <- file.path(res_dir,
                         res_files)
  names(res_paths) <- names(res_files)
  as.list(res_paths)
}

##' List all paths to results files for a study.
##'
##' @title List paths for all result files
##' @param study Study name
##' @param repls Replicate numbers
##' @param opmods Operating model numbers
##' @param estmods Estimation model names
##' @param root_dir Root directory
##' @return List where each element has \code{rdata} and \code{indexcsv} paths
##' @author John K Best
##' @export
all_res_file_paths <- function(study, repls, opmods, estmods, root_dir = ".") {
  res_dirs <- file.path(study_dir(study, root_dir),
                        "results",
                        repl_dir(repls))
  opmods <- paste0(study_file_base(study),
                   stringr::str_pad(opmods, 2, pad = 0))
  ## Names of files that should result from fits within each replicate
  res_files <- purrr::cross2(estmods, opmods) %>%
    purrr::map_chr(~ paste(.[2], .[1], sep = "_"))
  res_paths <- purrr::cross2(res_files,
                             res_dirs) %>%
    purrr::map_chr(~ file.path(.[2], .[1]))
  list(rdata = paste0(res_paths, ".Rdata"),
       indexcsv = paste0(res_paths, "_index.csv"))
}
