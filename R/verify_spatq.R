##' @title Get `dim` if available, `length` otherwise
##' @param x Object to be `dim`d
##' @return Vector (possible length 1) of the dimensions of `x`.
dim_or_len <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    ## Always return two dimensions
    d <- c(length(x), 1)
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
  verify_ident_fixef(data)
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
  datanames <- c("catch_obs", "area_swept",
                 "X_n", "X_w", "IX_n", "IX_w",
                 "Z_n", "Z_w", "IZ_n", "IZ_w",
                 "A_spat", "A_sptemp", "IA_spat", "IA_sptemp", "Ih",
                 "R_n", "R_w", "V_n", "V_w",
                 "A_qspat", "A_qsptemp", "spde")
  parnames <- c("beta_n", "beta_w", "gamma_n", "gamma_w",
                "omega_n", "omega_w", "epsilon1_n", "epsilon1_w",
                "lambda_n", "lambda_w", "eta_n", "eta_w",
                "phi_n", "phi_w", "psi1_n", "psi1_w",
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
    ddims$catch_obs[1] == ddims$area_swept[1]
    ddims$catch_obs[1] == ddims$X_n[1]
    ddims$catch_obs[1] == ddims$X_w[1]
    ddims$catch_obs[1] == ddims$X_n[1]
    ddims$catch_obs[1] == ddims$X_w[1]
    ddims$catch_obs[1] == ddims$A_spat[1]
    ddims$catch_obs[1] == ddims$A_sptemp[1]
    ddims$catch_obs[1] == ddims$R_n[1]
    ddims$catch_obs[1] == ddims$R_w[1]
    ddims$catch_obs[1] == ddims$V_n[1]
    ddims$catch_obs[1] == ddims$V_w[1]
    ddims$catch_obs[1] == ddims$A_qspat[1]
    ddims$catch_obs[1] == ddims$A_qsptemp[1]
  })

  ## Check that index components are the correct length
  stopifnot(exprs = {
    ddims$IX_n[1] == ddims$IX_w[1]
    ddims$IX_n[1] == ddims$IZ_n[1]
    ddims$IX_n[1] == ddims$IZ_w[1]
    ddims$IX_n[1] == ddims$IA_spat[1]
    ddims$IX_n[1] == ddims$IA_sptemp[1]
  })

  ## Check that matrices conform
  stopifnot(exprs = {
    ddims$X_n[2] == pdims$beta_n[1]
    ddims$X_w[2] == pdims$beta_w[1]
    ddims$Z_n[2] == pdims$gamma_n[1]
    ddims$Z_w[2] == pdims$gamma_w[1]
    ddims$A_spat[2] == pdims$omega_n[1]
    ddims$A_spat[2] == pdims$omega_w[1]
    ## epsilon1_* have one fewer column than years for sum-to-zero constraint
    ddims$A_sptemp[2] == pdims$epsilon1_n[1] * (pdims$epsilon1_n[2] + 1)
    ddims$A_sptemp[2] == pdims$epsilon1_w[1] * (pdims$epsilon1_w[2] + 1)
    ddims$R_n[2] == pdims$lambda_n[1]
    ddims$R_w[2] == pdims$lambda_w[1]
    ddims$V_n[2] == pdims$eta_n[1]
    ddims$V_w[2] == pdims$eta_w[1]
    ddims$A_qspat[2] == pdims$phi_n[1]
    ddims$A_qspat[2] == pdims$phi_w[1]
    ## psi1_* have one fewer column than years for sum-to-zero constraint
    ddims$A_sptemp[2] == pdims$psi1_n[1] * (pdims$psi1_n[2] + 1)
    ddims$A_sptemp[2] == pdims$psi1_w[1] * (pdims$psi1_w[2] + 1)
  })

  ## Check that all spatiotemporal processes have same number of columns
  stopifnot(exprs = {
    pdims$epsilon1_w[2] == pdims$epsilon1_n[2]
    pdims$psi1_n[2] == pdims$epsilon1_n[2]
    pdims$psi1_w[2] == pdims$epsilon1_n[2]
  })

  ## Check that index design matrices have same number of cols as data design
  ## matrices
  stopifnot(exprs = {
    ddims$IX_n[2] == ddims$X_n[2]
    ddims$IX_w[2] == ddims$X_w[2]
    ddims$IX_n[2] == ddims$X_n[2]
    ddims$IX_w[2] == ddims$X_w[2]
    ddims$IA_spat[2] == ddims$A_spat[2]
    ddims$IA_sptemp[2] == ddims$A_sptemp[2]
  })

  ## Check that number of index locations is consistent for each component of
  ## linear predictors
  stopifnot(exprs = {
    ddims$Ih[1] == ddims$IX_n[1]
    ddims$Ih[1] == ddims$IX_w[1]
    ddims$Ih[1] == ddims$IZ_n[1]
    ddims$Ih[1] == ddims$IZ_w[1]
    ddims$Ih[1] == ddims$IA_spat[1]
    ddims$Ih[1] == ddims$IA_sptemp[1]
  })

  ## Check that random effects parameters are correct length
  stopifnot(exprs = {
    pdims$log_xi[1] == 4
    pdims$log_kappa[1] == 8
    pdims$log_tau[1] == 8
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
  ## Check that iid random effects parameters are map'd if corresponding fields
  ## are
  ranef_names <- c("gamma_n", "gamma_w", "eta_n", "eta_w")
  for (nm in names(map)) {
    if (nm %in% ranef_names) {
      idx <- which(ranef_names == nm)
      stopifnot(exprs = {
        ## Check that map includes the field parametes and the correct indices
        ## are NA
        !is.null(map$log_xi) && is.na(map$log_xi[idx])
      })
      ## if (any(parameters[[nm]] != 0)) {
      ##   warn(nm, " is map'd but not all zeros.")
      ## }
    }
  }

  ## Check that spat/sptemp field parameters are map'd if corresponding fields
  ## are
  spattemp_names <- c("omega_n", "omega_w", "epsilon1_n", "epsilon1_w",
                      "phi_n", "phi_w", "psi1_n", "psi1_w")
  for (nm in names(map)) {
    ## Only check spat-temp field parameters if spat-temp fields are map'd.
    ## Allows all parameters to be included in `map` even when not actually
    ## map'd.
    if ((nm %in% spattemp_names) && all(is.na(map[[nm]]))) {
      idx <- which(spattemp_names == nm)
      stopifnot(exprs = {
        ## Check that map includes the field parameters and the correct indices
        ## are NA
        !is.null(map$log_tau) && is.na(map$log_tau[idx])
        !is.null(map$log_kappa) && is.na(map$log_kappa[idx])
      })
      if (warn_not_zero && any(parameters[[nm]] != 0)) {
        warning(nm, " is map'd but not all zeros.")
      }
    }
  }
  return(TRUE)
}

##' Check that the combined fixed effect design matrices for abundance and
##' catchability are full rank.
##'
##' @title Verify that fixed effects are identifiable
##' @param data Data list, as produced by \code{prepare_data}
##' @return TRUE if it passes, error otherwise
##' @author John Best
##' @export
verify_ident_fixef <- function(data) {
  fix_n <- cbind(data$X_n, data$R_n)
  np_n <- ncol(fix_n)
  rank_n <- Matrix::rankMatrix(fix_n)
  if (rank_n < np_n) stop("Numbers density fixed effects are not identifiable")

  fix_w <- cbind(data$X_w, data$R_w)
  np_w <- ncol(fix_w)
  rank_w <- Matrix::rankMatrix(fix_w)
  if (rank_w < np_w) stop("Weight per group fixed effects are not identifiable")

  return(TRUE)
}
