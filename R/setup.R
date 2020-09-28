##' Convenience constructor when data follows simulation study directory
##' structure. See also \code{\link{new_spatqsetup}}.
##'
##' @title Construct a spatqsetup from a simulated data set
##' @param repl Replicate number
##' @param sc Scenario; "pref", "spat", or "combo"
##' @param sub_df Data frame indicating subsampling strategy; see
##'   \code{\link{subsample_catch}}
##' @param root_dir Directory to load data from
##' @param max_T Last year of data to include
##' @param index_step Step for the index grid
##' @param spec_estd List of logicals indicating which parameters are to be
##'   estimated, as output \code{\link{specify_estimated}}
##' @return A list with elements \code{data}, \code{parameters}, \code{map}, and
##'   \code{random}, ready to be used as the arguments to
##'   \code{\link{prepare_adfun}}.
##' @author John K Best
##' @export
spatqsetup_sim <- function(repl, sc, sub_df = NULL,
                           root_dir = ".", max_T = NULL,
                           index_step = 5,
                           spec_estd = specify_estimated()) {
  ## Read in data
  catch_df <- read_catch(repl, sc, root_dir)
  if (!is.null(max_T)) {
    catch_df <- dplyr::filter(catch_df, time <= max_T)
    attr(catch_df, "T") <- max_T
  }
  ## Subset observations
  catch_df <- subsample_catch(catch_df, sub_df)

  return(spatqsetup(catch_df, spec_estd, index_step))
}

##' Checks that \code{data}, \code{parameters}, and \code{map} are consistent
##' then puts them in a list for use by \code{\link{prepare_adfun}}
##'
##' @title Set up inputs for spatq model
##' @param data Data for spatq model, as produced by \code{\link{prepare_data}}
##' @param parameters Parameter list as from \code{\link{prepare_pars}}
##' @param map Map list as from \code{\link{prepare_map}}
##' @param random Random vector as from \code{\link{prepare_random}}
##' @return A \code{spatqsetup} object ready to be used to construct a spatq
##'   objective function
##' @author John K Best
##' @export
new_spatqsetup <- function(data, parameters, map, random) {
  verify_spatq(data, parameters, map)
  structure(list(data = data,
                 parameters = parameters,
                 map = map,
                 random = random),
            class = "spatqsetup")
}

##' @describeIn new_spatqsetup Convenient constructor for model setups
##' @param catch_df Data frame with catch observations, as from
##'   \code{\link{read_catch}}
##' @param spatqspec Model specification, \code{\link{specify_estimated}}
##' @param index_step Index grid step size \code{\link{create_index_df}}
##' @export
spatqsetup <- function(catch_df, spatqspec, index_step) {
  ## Get number of years represented in catch data
  T <- length(unique(catch_df$time))

  ## Create index integration reference
  index_df <- create_index_df(step = index_step, T = T)

  ## Discretize space
  mesh <- generate_mesh()
  fem <- generate_fem(mesh)

  ## Prepare model specification components
  data <- prepare_data(catch_df, index_df, mesh, fem)
  parameters <- prepare_pars(data, mesh)
  map <- prepare_map(parameters,
                     spec = spatqspec)
  random <- prepare_random(map)

  return(new_spatqsetup(data, parameters, map, random))
}

##' Not intended to be used when the underlying data change.
##'
##' @title Update a model setup
##' @param setup Existing setup object
##' @param obj Objective function after optimizing
##' @param newspec Next estimation specification
##' @return Updated \code{\link{spatqsetup}}
##' @author John K Best
##' @export
update_setup <- function(setup, obj, newspec) {
  setup <- update_parameters(setup, obj)
  setup <- update_map(setup, newspec)
  setup <- update_random(setup)
  return(setup)
}

##' @describeIn update_setup Update \code{parameters}
##' @export
update_parameters <- function(setup, obj) {
  ## Extract the best values from the previous fit
  newpars <- gather_nvec(obj$env$last.par.best)
  ## Only need to iterate over parameters that were updated
  for (nm in names(newpars)) {
    ## Need to deal with parameters that are partially map'd
    if (nm %in% names(setup$map)) {
      ## Which parameter elements were not estimated?
      map_idx <- !is.na(setup$map[[nm]])
      ## Don't replace if there aren't the right number of new parameter values
      if (sum(map_idx) != length(newpars[[nm]])) {
        stop("New parameter values wrong length wrt map")
      }
      setup$parameters[[nm]][map_idx] <- newpars[[nm]]
    } else {
      ## If the parameter doesn't appear in `map`, replace the entire vector
      setup$parameters[[nm]] <- newpars[[nm]]
    }
  }

  ## Return entire modified `setup` object
  return(setup)
}

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
