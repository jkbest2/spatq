##' Save the results of a fit.
##'
##' @title Save fit results
##' @param studyspec A \code{spatq_studyspec} specification
##' @param fit A \code{\link{spatq_fit}} fit result
##' @param lpb The contents of \code{obj$env$last.par.best}, gathered into a
##'   named list
##' @param rep A parameter report from \code{\link{report_spatq}}
##' @param sdr An SD Report from \code{\link{sdreport_spatq}}
##' @param root_dir Root directory
##' @return The name of the saved file
##' @author John K Best
##' @export
save_fit <- function(studyspec, fit, lpb, rep, sdr, root_dir = ".") {
  rdata_path <- res_file_paths(studyspec$study,
                               studyspec$repl,
                               studyspec$opmod,
                               studyspec$estmod,
                               root_dir)
  saveRDS(list(spec = studyspec,
               fit = fit,
               lpb = lpb,
               rep = rep,
               sdr = sdr),
          rdata_path$rdata)
  invisible(rdata_path$rdata)
}

##' Save the index results, scaled appropriately. Also includes error estimates.
##'
##' @title Save true and estimated indices to a CSV
##' @param studyspec Simulation study specification
##' @param sdr  An SD report from \code{\link{sdreport_spatq}}
##' @param max_T End year
##' @param feather Save as feather or CSV?
##' @return The CSV file name
##' @author John K Best
##' @export
save_index <- function(studyspec, sdr, max_T = 15, feather = TRUE) {
  ## Read true population state and calculate index
  true_index <- read_popstate(study = studyspec$study,
                              repl = studyspec$repl,
                              opmod = studyspec$opmod,
                              root_dir = studyspec$root_dir) %>%
    dplyr::rename(raw_true = pop) %>%
    dplyr::filter(year <= max_T) %>%
    dplyr::mutate(index_true = rescale_index(raw_true)$index)

  if (!("fail" %in% names(sdr))) {
    ## Organize details for estimated index
    which_index <- which(names(sdr$value) == "Index")
    est_index <- tibble::tibble(study = studyspec$study,
                                repl = studyspec$repl,
                                opmod = studyspec$opmod,
                                estmod = studyspec$estmod,
                                year = 1:max_T,
                                raw_est = sdr$value[which_index],
                                index_est = rescale_index(raw_est)$index,
                                raw_unb = sdr$unbiased$value[which_index],
                                index_unb = rescale_index(raw_unb)$index,
                                raw_sd = sdr$sd[which_index],
                                index_sd = rescale_index(raw_est, raw_sd)$sd,
                                raw_unb_sd = sdr$unbiased$sd,
                                unb_sd = rescale_index(raw_unb, raw_unb_sd)$sd)
  } else {
    est_index <- tibble::tibble(study = studyspec$study,
                                repl = studyspec$repl,
                                opmod = studyspec$opmod,
                                estmod = studyspec$estmod,
                                year = 1:max_T,
                                raw_est = rep(NA, max_T),
                                index_est = rep(NA, max_T),
                                raw_unb = rep(NA, max_T),
                                index_unb = rep(NA, max_T),
                                raw_sd = rep(NA, max_T),
                                index_sd = rep(NA, max_T),
                                raw_unb_sd = rep(NA, max_T),
                                unb_sd = rep(NA, max_T))

  }

  ## Join and write
  index_df <- dplyr::left_join(est_index, true_index, by = "year")
  index_path <- res_file_paths(studyspec$study,
                               studyspec$repl,
                               studyspec$opmod,
                               studyspec$estmod,
                               studyspec$root_dir)
  if (feather) {
    flnm <- index_path$index_feather
    arrow::write_feather(index_df, flnm)
  } else {
    flnm <- index_path$index_csv
    readr::write_csv(index_df, flnm)
  }
  flnm
}
