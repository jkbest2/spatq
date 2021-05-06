##' Save the results of a fit.
##'
##' @title Save fit results
##' @param spec A \code{spatq_spec} specification
##' @param fit A \code{\link{spatq_fit}} fit result
##' @param lpb The contents of \code{obj$env$last.par.best}, gathered into a
##'   named list
##' @param rep A parameter report from \code{\link{report_spatq}}
##' @param sdr An SD Report from \code{\link{sdreport_spatq}}
##' @return The name of the saved file
##' @author John K Best
save_fit <- function(spec, fit, lpb, rep, sdr) {
  saveRDS(list(spec = spec,
               fit = fit,
               lpb = lpb,
               rep = rep,
               sdr = sdr),
          spec$Rdata)
  invisible(spec$Rdata)
}

##' Save the index results, scaled appropriately. Also includes error estimates.
##'
##' @title Save true and estimated indices to a CSV
##' @param spec A \code{spatq_spec}
##' @param sdr  An SD report from \code{\link{sdreport_spatq}}
##' @param max_T End year
##' @return The CSV file name
##' @author John K Best
save_index <- function(spec, sdr, max_T = 15) {
  ## Read true population state and calculate index
  true_index <- read_popstate(study = spec$study,
                              repl = spec$repl,
                              opmod = spec$opmod,
                              root_dir = ".") %>%
    dplyr::rename(year = time,
                  raw_true = pop) %>%
    dplyr::filter(year <= max_T) %>%
    dplyr::mutate(index_true = rescale_index(raw_true)$index,
                  study = spec$study)

  if (!("fail" %in% names(sdr))) {
    ## Organize details for estimated index
    which_index <- which(names(sdr$value) == "Index")
    est_index <- tibble::tibble(study = spec$study,
                                repl = spec$repl,
                                opmod = spec$opmod,
                                estmod = spec$estmod,
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
    est_index <- tibble::tibble(study = spec$study,
                                repl = spec$repl,
                                opmod = spec$opmod,
                                estmod = spec$estmod,
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

  ## Join and write to CSV file
  index_df <- dplyr::left_join(est_index, true_index, by = "year")
  readr::write_csv(index_df, spec$index)
  invisible(spec$index)
}
