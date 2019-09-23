##' @title Get `dim` if available, `length` otherwise
##' @param x Object to be `dim`d
##' @return Vector (possible length 1) of the dimensions of `x`.
dim_or_len <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    d <- length(x)
  }
  d
}

##' Check data, parameter, and map names, dimensions, and `map` specification.
##' Convenient wrapper around `verify_spatq_names`, `verify_spatq_dims`, and
##' `verify_spatq_map`.
##'
##' @title Check `spatq` inputs
##' @param data Data list
##' @param parameters Initial parameter value list
##' @param map Parameter map list
##' @return `TRUE` if all checks pass
##' @export
verify_spatq <- function(data, parameters, map = list()) {
  verify_spatq_names(data, parameters, map)
  verify_spatq_dims(data, parameters, map)
  verify_spatq_map(parameters, map)
}

##' Check that all the required data and parameter names are included. Currently
##' doesn't check `map` names.
##'
##' @title Check all `spatq` components are present
##' @param data Data list for `MakeADFun`
##' @param parameters Initial parameter value list
##' @param map Parameter map list
##' @return `TRUE` if all components present
##' @export
verify_spatq_names <- function(data, parameters, map = list()) {
  ## Check that all required elements are present
  datanames <- c("catch_obs", "X_n", "X_w", "A_spat", "A_sptemp",
                 "R_n", "R_w", "A_qspat", "A_qsptemp", "spde")
  parnames <- c("beta_n", "beta_w", "omega_n", "omega_w",
                "epsilon_n", "epsilon_w",
                "lambda_n", "lambda_w", "phi_n", "phi_w",
                "psi_n", "psi_w",
                "log_kappa", "log_tau", "log_sigma")
  ## print(parameters)
  stopifnot(exprs = {
    all(datanames %in% names(data))
    all(parnames %in% names(parameters))
  })
  return(TRUE)
}

##' Check that all components are of correct dimensions and conform.
##'
##' @title Check `spatq` component dimensions
##' @param data Data list for `MakeADFun`
##' @param parameters Initial parameter value list
##' @param map Parameter map list
##' @return `TRUE` if all dimension checks pass
##' @export
verify_spatq_dims <- function(data, parameters, map = list()) {
  ddims <- lapply(data, dim_or_len)
  pdims <- lapply(parameters, dim_or_len)
  mdims <- lapply(map, dim_or_len)

  ## Check that all vector multiplies will result in vector same length as
  ## number of obserations
  stopifnot(exprs = {
    ddims$catch_obs == ddims$X_n[1]
    ddims$catch_obs == ddims$X_w[1]
    ddims$catch_obs == ddims$A_spat[1]
    ddims$catch_obs == ddims$A_sptemp[1]
    ddims$catch_obs == ddims$R_n[1]
    ddims$catch_obs == ddims$R_w[1]
    ddims$catch_obs == ddims$A_qspat[1]
    ddims$catch_obs == ddims$A_qsptemp[1]
  })

  ## Check that matrices conform
  stopifnot(exprs = {
    ddims$X_n[2] == pdims$beta_n
    ddims$X_w[2] == pdims$beta_w
    ddims$A_spat[2] == pdims$omega_n
    ddims$A_spat[2] == pdims$omega_w
    ddims$A_sptemp[2] == prod(pdims$epsilon_n)
    ddims$A_sptemp[2] == prod(pdims$epsilon_w)
    ddims$R_n[2] == pdims$lambda_n
    ddims$R_w[2] == pdims$lambda_w
    ddims$A_qspat[2] == pdims$phi_n
    ddims$A_qspat[2] == pdims$phi_w
    ddims$A_qsptemp[2] == prod(pdims$psi_n)
    ddims$A_qsptemp[2] == prod(pdims$psi_w)
  })

  ## Check that all spatiotemporal processes have same number of years
  stopifnot(exprs = {
    pdims$epsilon_w[2] == pdims$epsilon_n[2]
    pdims$psi_n[2] == pdims$epsilon_n[2]
    pdims$psi_w[2] == pdims$epsilon_n[2]
  })

  ## Check that log_kappa and log_tau are each length 8
  stopifnot(exprs = {
    pdims$log_kappa == 8
    pdims$log_tau == 8
  })

  lapply(names(map), function(mn) {
    stopifnot(length(map[[mn]]) == length(parameters[[mn]]))
  })

  return(TRUE)
}

##' Check that spatial or spatiotemporal field parameters are map'd off if the
##' corresponding field coefficients are map'd off.
##'
##' @title Check that `map` is self-consistent
##' @param parameters Initial parameter value list
##' @param map Parameter map list
##' @param warn_not_zero Raise a warning if spatial or spatiotemporal components
##'   are map'd off but not all equal to zero.
##' @return `TRUE` if all checks pass.
##' @export
verify_spatq_map <- function(parameters, map, warn_not_zero = TRUE) {
  ## Check that spat/sptemp field parameters are map'd if corresponding fields
  ## are
  spattemp_names <- c("omega_n", "omega_w", "epsilon_n", "epsilon_w",
                      "phi_n", "phi_w", "psi_n", "psi_w")
  for (nm in names(map)) {
    if (nm %in% spattemp_names) {
      idx <- which(spattemp_names == nm)
      stopifnot(exprs = {
        ## Check that map includes the field parametes and the correct indices
        ## are NA
        !is.null(map$log_tau) && is.na(map$log_tau[idx])
        !is.null(map$log_kappa) && is.na(map$log_kappa[idx])
      })
      if (any(parameters[[nm]] != 0)) {
        warn(nm, " is map'd but not all zeros.")
      }
    }
  }
  return(TRUE)
}