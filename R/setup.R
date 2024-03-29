##' Convenience constructor when data follows simulation study directory
##' structure. See also \code{\link{new_spatq_setup}}.
##'
##' @title Construct a spatq_setup from a simulated data set
##' @param study Simulation study name
##' @param repl Replicate number
##' @param opmod Operating model
##' @param estmod Estimation model
##' @param sub_df Data frame indicating subsampling strategy; see
##'   \code{\link{subsample_catch}}
##' @param root_dir Directory to load data from
##' @param max_T Last year of data to include
##' @param index_step Step for the index grid
##' @param spec_estd List of logicals indicating which parameters are to be
##'   estimated, as output \code{\link{specify_estimated}}
##' @param mesh_resolution Mesh resolution to pass to
##'   \code{\link{generate_mesh}}. Must be "coarse", "medium", or "fine".
##' @param init_fixef Use parameter estimates from model with only fixed effects
##'   as initial values?
##' @param ... Additional arguments to pass to \code{\link{init_fixef}}, which
##'   are then passed to \code{\link{spatq_fit}} for the fixed effect model fit
##' @return A list with elements \code{data}, \code{parameters}, \code{map}, and
##'   \code{random}, ready to be used as the arguments to
##'   \code{\link{prepare_adfun}} or passed to \code{\link{spatq_obj}}
##' @author John K Best
##' @export
spatq_simsetup <- function(study, repl, opmod, estmod = NULL,
                           sub_df = NULL, spec_estd = NULL,
                           root_dir = ".", max_T = NULL,
                           index_step = 5,
                           mesh_resolution = "coarse",
                           init_fixef = TRUE, ...) {
  if (is.null(sub_df)) sub_df <- em_subsample(estmod)
  if (is.null(spec_estd)) spec_estd <- em_estd(estmod)

  ## Read in data
  catch_df <- read_catch(study, repl, opmod, root_dir)
  if (!is.null(max_T)) {
    catch_df <- dplyr::filter(catch_df, year <= max_T)
    attr(catch_df, "T") <- max_T
  }
  ## Subset observations
  catch_df <- subsample_catch(catch_df, sub_df)

  ## For "design*" EMs, just return the data frame
  if (!is.null(estmod) && estmod %in% c("design", "design_all")) {
    class(catch_df) <- c("spatq_designsetup", class(catch_df))
    return(catch_df)
  }

  setup <- spatq_setup(catch_df, spec_estd, index_step, init_fixef, ...)
  attr(setup, "repl") <- repl
  attr(setup, "opmod") <- opmod
  attr(setup, "estmod") <- estmod
  class(setup) <- c("spatq_simsetup", "spatq_setup")
  return(setup)
}

##' Checks that \code{data}, \code{parameters}, and \code{map} are consistent
##' then puts them in a list for use by \code{\link{prepare_adfun}}
##'
##' @title Set up inputs for spatq model
##' @param data Data for spatq model, as produced by \code{\link{prepare_data}}
##' @param parameters Parameter list as from \code{\link{prepare_pars}}
##' @param map Map list as from \code{\link{prepare_map}}
##' @param random Random vector as from \code{\link{prepare_random}}
##' @return A \code{spatq_setup} object ready to be used to construct a spatq
##'   objective function
##' @author John K Best
##' @export
new_spatq_setup <- function(data, parameters, map, random) {
  verify_spatq(data, parameters, map)
  structure(list(
    data = data,
    parameters = parameters,
    map = map,
    random = random
  ),
  class = "spatq_setup"
  )
}

##' @describeIn new_spatq_setup Convenient constructor for model setups
##' @param catch_df Data frame with catch observations, as from
##'   \code{\link{read_catch}}
##' @param spatq_spec Model specification from \code{\link{specify_estimated}}
##' @param index_step Index grid step size \code{\link{create_index_df}}
##' @param init_fixef Use parameter estimates from model with only fixed effects
##'   as initial values?
##' @param mesh_resolution Mesh resolution to pass to
##'   \code{\link{generate_mesh}}. Must be "coarse", "medium", or "fine".
##' @param ... Additional arguments to pass to \code{\link{init_fixef}}, which
##'   are then passed to \code{\link{spatq_fit}} for the fixed effect model fit
##' @export
spatq_setup <- function(catch_df, spatq_spec, index_step, init_fixef = TRUE, mesh_resolution = "coarse", ...) {
  ## Get number of years represented in catch data
  T <- length(unique(catch_df$year))

  ## Create index integration reference
  index_df <- create_index_df(step = index_step, T = T)

  ## Discretize space
  mesh <- generate_mesh(resolution = mesh_resolution)
  fem <- generate_aniso_fem(mesh)

  ## Prepare model specification components
  data <- prepare_data(catch_df, index_df, mesh, fem, spatq_spec$obs_lik)
  parameters <- prepare_pars(data, mesh)
  map <- prepare_map(parameters,
    spec = spatq_spec
  )
  random <- prepare_random(map)

  setup <- new_spatq_setup(data, parameters, map, random)

  if (init_fixef) {
    setup <- init_fixef(setup, spatq_spec, ...)
    verify_spatq(setup$data, setup$parameters, setup$map)
  }

  return(setup)
}

##' Fit the specified model without any random effects and use the resulting
##' parameter estimates as initial values.
##'
##' @title Initial estimates of fixed effect parameters
##' @param setup A \code{\link{spatq_setup}}
##' @param spec A model specification, as from \code{\link{specify_estimated}}
##' @param ... Further arguments to pass to \code{\link{spatq_fit}}
##' @return A \code{\link{spatq_setup}} with parameters of fixed-effects model
##'   initialized to their MLEs
##' @author John Best
##' @export
init_fixef <- function(setup, spec, ...) {
  fix_spec <- specify_estimated(
    beta = TRUE,
    lambda = !attr(setup$parameters, "map_lambda"),
    obs_lik = setup$data$obs_lik
  )
  fix_setup <- update_setup(setup, setup$parameters, fix_spec)
  fix_obj <- spatq_obj(fix_setup, runSymbolicAnalysis = FALSE, silent = TRUE)
  fix_fit <- spatq_fit(fix_obj, ...)
  setup <- update_setup(setup, fix_fit, spec)
}

##' Not intended to be used when the underlying data change.
##'
##' @title Update a model setup
##' @param setup Existing setup object
##' @param newparobj Object containing updated parameter values, to be extracted
##'   using \code{\link{get_newpars}}
##' @param newspec Next estimation specification
##' @return Updated \code{\link{spatq_setup}}
##' @author John K Best
##' @export
update_setup <- function(setup, newparobj, newspec) {
  setup <- update_parameters(setup, newparobj)
  setup <- update_map(setup, newspec)
  setup <- update_random(setup)
  return(setup)
}

##' @describeIn update_setup Update \code{parameters}
##' @export
update_parameters <- function(setup, newparobj) {
  ## Extract the best values from the previous fit as a named list
  newpars <- get_newpars(newparobj)
  ## Update each parameter vector using `newpars`. Use a `for` loop because need
  ## to iterate over `names(setup$parameters)`, so `lapply` returns an unnamed
  ## list.
  for (nm in names(setup$parameters)) {
    setup$parameters[[nm]] <- update_onepar(
      setup$parameters[[nm]],
      newpars[[nm]],
      setup$map[[nm]]
    )
  }

  ## Return entire modified `setup` object
  return(setup)
}

##' Extract updated parameter values from a previous fit, provided as a named
##' vector, list, \code{spatq_fit}, or \code{spatq_obj} that has already been
##' optimized.
##'
##' @title Get updated parameter values
##' @param obj object to extract new parameter values from
##' @return named list with parameter values to update
##' @author John K Best
##' @export
get_newpars <- function(obj) UseMethod("get_newpars")
##' @method get_newpars list
##' @export
get_newpars.list <- function(obj) obj
##' @method get_newpars spatq_fit
##' @export
get_newpars.spatq_fit <- function(obj) gather_nvec(obj$par)
##' @method get_newpars spatq_obj
##' @export
get_newpars.spatq_obj <- function(obj) gather_nvec(obj$env$last.par.best)
##' @method get_newpars sdreport
##' @export
get_newpars.sdreport <- function(obj) gather_nvec(c(obj$par.fixed, obj$par.random))
##' @method get_newpars numeric
## Vectors inherit "numeric"
##' @export
get_newpars.numeric <- function(obj) gather_nvec(obj)
## Handle the NULL case
##' @method get_newpars NULL
##' @export
get_newpars.NULL <- function(obj) NULL

##' @describeIn update_setup Update \code{map}
##' @export
update_map <- function(setup, newspec) {
  setup$map <- prepare_map(setup$parameters, newspec)
  return(setup)
}

##' @describeIn update_setup Update \code{random}
##' @export
update_random <- function(setup) {
  setup$random <- prepare_random(setup$map)
  setup
}

##' Replace parameter values with better ones, if available.
##'
##' @title Update a single parameter vector
##' @param currvals Vector of current parameter values
##' @param newvals Values to update
##' @param mapvec Map, indexing the new values into their places in the
##'   parameter vector, with NAs for elements that were not estiamted
##' @return The updated parameter vector
##' @author John K Best
##' @importFrom utils tail
##' @export
update_onepar <- function(currvals, newvals = NULL, mapvec = NULL) {
  if (is.null(newvals)) {
    ## If no new values, use old values
    vals <- currvals
  } else if (is.null(mapvec)) {
    ## If not map'd, just replace the vector
    vals <- newvals
  } else if (all(is.na(mapvec))) {
    ## If all map'd and new values provided, use new values
    vals <- newvals
  } else {
    ## Need to
    if (length(currvals) < length(mapvec)) {
      stop("Length of currvals must be greater than mapvec")
    }

    ## Going to use for indexing and need to check `max`
    mapvec <- as.integer(mapvec)

    if (max(mapvec, na.rm = TRUE) > length(newvals)) {
      ## If new parameters are to be estimated, use the current values
      mapvec[mapvec > length(newvals)] <- NA
    }

    ## If not estimated, use current value. Otherwise use value corresponding to
    ## the map index. Need this to account for parameters map'd to take same
    ## values (e.g. log_kappa map is c(1, 2, NA, NA, 1, 2, NA, NA)).
    vals <- ifelse(is.na(mapvec), currvals, newvals[mapvec])
  }

  ## Used when parameter dimension changes, e.g. adding years or a new fleet. In
  ## that case you won't have new parameter estimates, but we assume that
  ## they're added on the end of the parameter vector so we fill in the tail
  ## using the current (default) values.
  if (length(vals) < length(currvals)) {
    vals <- append(vals, tail(currvals, length(currvals) - length(vals)))
  }

  return(vals)
}
