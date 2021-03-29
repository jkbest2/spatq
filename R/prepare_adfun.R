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
  structure(obj, class = "spatq_obj")
}

##' Prepare an objective function using \code{\link[TMB]{MakeADFun}}.
##'
##' @title Contstruct a spatq objective function
##' @param setup a \code{\link{spatq_setup}} object
##' @param runSymbolicAnalysis Use Metis reorderings? (Requires special
##'   installation of TMB; see documentation.)
##' @param normalize Normalize GMRF likelihoods?
##' @param ... Additional arguments to pass to \code{\link[TMB]{MakeADFun}}
##' @return a spatq_obj object suitable for passing to \code{\link{fit_spatq}}
##' @author John K Best
##' @export
spatq_obj <- function(setup, runSymbolicAnalysis = TRUE, normalize = TRUE, ...) {
  setup$data$proc_switch <- prepare_proc_switch(setup$random)
  setup$data$norm_flag <- normalize
  verify_spatq(setup$data, setup$parameters, setup$map)
  obj <- TMB::MakeADFun(data = setup$data,
                        parameters = setup$parameters,
                        map = setup$map,
                        random = setup$random,
                        ...)
  if (!normalize)
    obj <- TMB::normalize(obj, flag = "incl_data", value = FALSE)
  if (runSymbolicAnalysis & length(setup$random) > 0)
    TMB::runSymbolicAnalysis(obj)
  structure(obj, class = "spatq_obj")
}

##' Read in a simulated data set and construct a TMB ADFun for model fitting.
##'
##' @title Create ADFun from simulated data
##' @param study Simulation study
##' @param repl Replicate number
##' @param opmod Operating model
##' @param sub_df Data frame indicating subsampling strategy; see
##'   \code{subsample_catch}
##' @param root_dir Directory to load data from
##' @param max_T Last year of data to include
##' @param index_step Index grid step size
##' @param spec_estd List of logicals indicating which parameters are to be
##'   estimated, as output \code{\link{specify_estimated}}
##' @param ... Additional options to pass to \code{prepare_adfun}
##' @return A TMB ADFun suitable for optimization
##' @author John Best
##' @export
make_sim_adfun <- function(study, repl, opmod, sub_df = NULL,
                           root_dir = ".", max_T = NULL,
                           index_step,
                           spec_estd = specify_estimated(), ...) {
  setup <- spatq_simsetup(study,
                          repl,
                          opmod,
                          sub_df,
                          root_dir,
                          max_T,
                          index_step,
                          spec_estd)
  prepare_adfun(setup$data, setup$parameters,
                setup$map, setup$random, ...)
}

