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
    counterpref = "counterpref",
    habq = "habq",
    bycatch = "bycatch",
    moverate = "moverate"
  )
  file.path(root_dir, sdir)
}

##' @describeIn study_dir The file base for each study
study_file_base <- function(study) {
  study_file <- switch(study,
    qdevscaling = "qdevscale_",
    sharedq = "sharedq_",
    prefintensity = "prefintensity_",
    densdepq = "densdepq_",
    counterpref = "counterpref_",
    habq = "habq_",
    bycatch = "bycatch_",
    moverate = "moverate_"
  )
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
  prep <- switch(study,
    qdevscaling = ,
    sharedq = ,
    prefintensity = ,
    densdepq = ,
    habq = ,
    bycatch = ,
    moverate = paste0(study_file, "prep.h5"),
    counterpref = "../../prep.h5"
  )
  opmod <- stringr::str_pad(opmod, 2, pad = "0")
  list(
    catch_csv = paste0(study_file, opmod, "_catch.csv"),
    catch_feather = paste0(study_file, opmod, "_catch.feather"),
    pop_csv = paste0(study_file, opmod, "_popstate.csv"),
    pop_feather = paste0(study_file, opmod, "_popstate.feather"),
    pop_h5 = paste0(study_file, opmod, "_popstate.h5"),
    prep = prep
  )
}

##' Get the paths to each simulation file. ' File names are
## ${study_dir}/repl_${repl}/${study_opmod}_catch.csv File ' names are
## repl_$repl/catch_$repl_$sc.csv, with $repl padded to two digits.
##' @title Paths to simulated data sets
##' @param study Study
##' @param repl Replicate number
##' @param opmod Operating model
##' @param root_dir Where to start
##' @return A list of file paths with elements \code{catch_csv},
##'   \code{catch_feather}, \code{pop_csv}, \code{pop_feather}, and \code{poph5}
##' @author John K Best
##' @export
sim_file_paths <- function(study, repl, opmod, root_dir = ".") {
  files <- sim_file_names(study, opmod)
  paths <- file.path(
    study_dir(study, root_dir),
    repl_dir(repl),
    files
  )
  names(paths) <- names(files)
  as.list(paths)
}

##' @describeIn res_file_paths Create a results directory if necessary
##' @export
create_res_dir <- function(study, repl = NULL, root_dir = ".") {
  res <- file.path(study_dir(study, root_dir), "results")
  if (!dir.exists(res)) {
    dir.create(res, recursive = TRUE)
  }
  if (!is.null(repl)) {
    res_repl <- file.path(
      study_dir(study, root_dir),
      "results",
      repl_dir(repl)
    )
    ## If ${study_dir}/resutls/${repl_dir} doesn't exist, create it
    vapply(
      res_repl,
      function(rdir) {
        v <- TRUE
        if (!dir.exists(rdir)) {
          v <- dir.create(rdir, recursive = TRUE)
        }
        v
      }, TRUE
    )
  }
  invisible(TRUE)
}

##' @describeIn res_file_paths File names for each study, operating, and
##'   estimation model combination
res_file_names <- function(study, opmod, estmod) {
  opmod <- stringr::str_pad(opmod, 2, pad = 0)
  res_files <- paste0(
    study_file_base(study),
    opmod, "_",
    estmod,
    c(".Rdata", "_index.csv", "_index.feather")
  )
  names(res_files) <- c("rdata", "index_csv", "index_feather")
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
  res_dir <- file.path(
    study_dir(study, root_dir),
    "results",
    repl_dir(repl)
  )
  res_files <- res_file_names(study, opmod, estmod)
  res_paths <- file.path(
    res_dir,
    res_files
  )
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
##' @return List where each element has \code{rdata}, \code{index_csv}, and
##'   \code{index_feather} paths
##' @author John K Best
##' @export
all_res_file_paths <- function(study, repls, opmods, estmods, root_dir = ".") {
  res_dirs <- file.path(
    study_dir(study, root_dir),
    "results",
    repl_dir(repls)
  )
  opmods <- paste0(
    study_file_base(study),
    stringr::str_pad(opmods, 2, pad = 0)
  )
  ## Names of files that should result from fits within each replicate
  res_files <- purrr::cross2(estmods, opmods) %>%
    purrr::map_chr(~ paste(.[2], .[1], sep = "_"))
  res_paths <- purrr::cross2(
    res_files,
    res_dirs
  ) %>%
    purrr::map_chr(~ file.path(.[2], .[1]))
  list(
    rdata = paste0(res_paths, ".Rdata"),
    index_csv = paste0(res_paths, "_index.csv"),
    index_feather = paste0(res_paths, "_index.feather")
  )
}

##' @title Construct the relative path to the index file output by a fit.
##' @param spec A \code{\link{spatq_simstudyspec}} object
##' @param filetype Index file type, either "feather"  or "csv"
##' @return The relative path to the index data file
##' @author John K Best
##' @export
index_path <- function(spec, filetype) UseMethod("index_path")
##' @export
index_path.spatq_simstudyspec <- function(spec, filetype) {
  ## Check filetype
  filetype %in% c("feather", "csv") || stop("filetype must be \"feather\" or \"csv\"")
  opmod <- stringr::str_pad(spec$opmod, 2, pad = 0)
  index_file <- paste0(
    study_file_base(spec$study),
    opmod, "_",
    spec$estmod, "_index.",
    filetype
  )

  file.path(
    study_dir(spec$study, spec$root_dir),
    "results",
    repl_dir(spec$repl),
    index_file
  )
}
##' @export
index_path.list <- function(spec, filetype) {
  spec <- spatq_simstudyspec(spec)
  index_path(spec, filetype)
}

##' @title Construct the relative path to the Rdata fit file
##' @param spec A \code{\link{spatq_simstudyspec}} object or list with its
##'   elements
##' @return The relative path to the Rdata file
##' @author John K Best
##' @export
rdata_path <- function(spec) UseMethod("rdata_path")
##' @export
rdata_path.spatq_simstudyspec <- function(spec) {
  opmod <- stringr::str_pad(spec$opmod, 2, pad = 0)
  rdata_file <- paste0(
    study_file_base(spec$study),
    opmod, "_",
    spec$estmod, ".Rdata"
  )

  file.path(
    study_dir(spec$study, spec$root_dir),
    "results",
    repl_dir(spec$repl),
    rdata_file
  )
}
##' @export
rdata_path.list <- function(spec) {
  spec <- spatq_simstudyspec(spec)
  rdata_path(spec)
}
