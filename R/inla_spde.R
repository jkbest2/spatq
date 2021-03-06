##' Assumes the simulation \eqn{[0, 100]×[0, 100]} domain.
##'
##' @title Define the domain boundary
##' @return An \code{inla.mesh.segment} giving the boundary of the 100 by 100
##'   domain
##' @author John Best
##' @export
domain_boundary <- function() {
  boundary <- INLA::inla.mesh.segment(matrix(c(0, 0,   100, 100, 0,
                                               0, 100, 100, 0,   0),
                                             ncol = 2))
}

##' Generate the coordinates of a regular grid on a \eqn{[0, 100]×[0, 100]}
##' domain with steps in both direction of \code{step}.
##'
##' @title Generate a location grid
##' @param step Distance between subsequent locations; defaults to 1
##' @return A two-column matrix with of locations in two dimensions, each
##'   \code{step} apart.
##' @author John Best
##' @export
loc_grid <- function(step = 1.0) {
  grid_start <- step / 2
  grid_end <- 100 - step / 2
  grid_seq <- seq(grid_start, grid_end, step)
  as.matrix(expand.grid(s1 = grid_seq,
                        s2 = grid_seq))
}

##' Generate the standard mesh used for simulation fits
##'
##' @title Generate INLA mesh
##' @return A \code{inla.mesh} object that can be passed to \code{generate_fem}.
##' @author John Best
##' @export
generate_mesh <- function() {
  ## Discretize spatial domain into mesh
  boundary <- domain_boundary()
  loc <- loc_grid(2.0)
  INLA::inla.mesh.2d(loc,
                     boundary = boundary,
                     offset = c(0.0, 30.0),
                     # Shortest correlation range should be ~30
                     max.edge = c(5, 20),
                     max.n = c(400, 100),
                     min.angle = c(30, 21),
                     cutoff = 5)
}

##' Generate the FEM matrices with appropriate names to pass to TMB as a
##' \code{spde_t}.
##'
##' @title Generate FEM matrices
##' @param mesh INLA mesh to generate FEM matrices from
##' @return List of sparse matrices suitable for passing to TMB.
##' @author John Best
##' @export
generate_fem <- function(mesh) {
  ## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended
  ## by Finn Lindgren
  fem <- INLA::inla.mesh.fem(mesh)
  ## Rename elements of `fem` so that they are recognized as a `spde` object in
  ## TMB
  names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")
  fem
}

##' Generate a projection matrix to locations in \code{data_df}, possible
##' zeroing out some observations. Also group by year for a spatiotemporal
##' effect.
##'
##' @title Generate a projection matrix
##' @param mesh INLA mesh to project
##' @param data_df Data frame with \code{s1} and \code{s2} columns as
##'   coordinates to project to
##' @param vessel_idx Vessel index to project to; others receive zero weight and
##'   are dropped
##' @param group Year of observation for spatiotemporal effect
##' @param zero Is this an empty projection matrix? (May be useful if e.g. the
##'   random effects parameters associated with it are map'd to zeros.)
##' @return A sparse projection matrix
##' @author John Best
##' @export
generate_projection <- function(mesh, data_df, vessel_idx = NULL, group = NULL,
                                zero = FALSE) {
  ## If effect is map'd, return an all-zero projection matrix. Should save some
  ## unnecessary multiplications?
  if (zero) return(generate_empty_projection(mesh, data_df, group))
  if (is.null(vessel_idx)) {
    wts <- NULL
  } else {
    wts <- data_df$vessel_idx %in% vessel_idx
  }
  loc <- as.matrix(data_df[, c("s1", "s2")])
  A <- INLA::inla.spde.make.A(mesh, loc, group = group, weights = wts)
  Matrix::drop0(A)
}

##' When spatial or spatiotemporal effects are map'd to zero, pass an empty
##' projection matrix. Not exported, use `zero = TRUE` in `generate_projection`
##' instead.
##'
##' @title Generate an empty project matrix
##' @param mesh INLA mesh to project
##' @param data_df Data frame with \code{s1} and \code{s2} columns as
##'   coordinates to project to
##' @param group Year of observation for spatiotemporal effect
##' @return Sparse \code{Matrix::Matrx} of the appropriate dimensions, but
##'   filled with zeros
##' @author John Best
generate_empty_projection <- function(mesh, data_df, group = NULL) {
  n_obs <- nrow(data_df)
  if (!is.null(group)) {
    n_years <- length(unique(group))
  } else {
    n_years <- 1
  }
  Matrix::Matrix(0, nrow = n_obs, ncol = mesh$n * n_years)
}

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. The kappa parameter is \code{sqrt(8) / rho}.
##'
##' @describeIn pars_tau Calculate kappa from rho
##' @export
pars_kappa <- function(rho) {
  sqrt(8) / rho
}

##' The SPDE construction uses parameters kappa and tau, while the usual Matern
##' parameters are sigma and rho. These functions allow converting between the
##' parameterizations.
##'
##' The tau parameter rescales the precision matrix to a given marginal standard
##' deviation, while \code{kappa} controls the correlation decay rate. In particular,
##'
##' - \code{tau = 2 * sqrt(pi) * sig * kappa},
##' - \code{kappa = sqrt(8) / rho},
##' - \code{sig = tau / (2 * sqrt(pi) * kappa)}, and
##' - \code{rho = sqrt(8) / kappa}.
##'
##' All assume that the Matern smoothness parameter is 1 and the domain is
##' two-dimensional
##'
##' @title Calculate tau from rho and sigma^2
##' @param sig Marginal standard deviation
##' @param rho Correlation range parameter
##' @param kappa Correlation decay parameter
##' @param tau Ratio of precision and correlation range
##' @return Parameter value
##'
##' @author John Best
##' @export
pars_tau <- function(sig, rho) {
  2 * sig * sqrt(pi) * pars_kappa(rho)
}

##' @describeIn pars_tau Calculate the correlation range
##' @export
pars_rho <- function(kappa) {
  sqrt(8) / kappa
}

##' @describeIn pars_tau Calculate the marginal standard deviation
##' @export
pars_sig <- function(tau, kappa) {
  tau / (2 * sqrt(pi) * kappa)
}
