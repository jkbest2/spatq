##' Assumes the simulation \eqn{[0, 100]×[0, 100]} domain.
##'
##' @title Define the domain boundary
##' @return An \code{inla.mesh.segment} giving the boundary of the 100 by 100
##'   domain
##' @author John Best
##' @export
domain_boundary <- function() {
  boundary <- INLA::inla.mesh.segment(matrix(c(
    0, 0, 100, 100, 0,
    0, 100, 100, 0, 0
  ),
  ncol = 2
  ))
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
  as.matrix(expand.grid(
    s1 = grid_seq,
    s2 = grid_seq
  ))
}

##' Generate the standard mesh used for simulation fits
##'
##' Provides three options for standard meshes. All include the same constraints
##' on max edge length and minimum vertex angle, as well as offset. The "coarse"
##' mesh uses 404 vertices with a \code{cutoff} of 5, the "medium" mesh uses
##' 1437 vertices with a \code{cutoff} of 2, and the "fine" mesh uses 3038
##' vertices with a \code{cutoff} of 1.
##' @title Generate INLA mesh
##' @param resolution One of "coarse", "medium", or "fine"
##' @return A \code{inla.mesh} object that can be passed to \code{generate_fem}.
##' @author John Best
##' @export
generate_mesh <- function(resolution = "coarse") {
  ## Discretize spatial domain into mesh
  boundary <- domain_boundary()
  loc <- loc_grid(2.0)
  if (resolution == "coarse") {
    mesh <- INLA::inla.mesh.2d(loc,
      boundary = boundary,
      offset = c(0.0, 30.0),
      ## Shortest correlation range should be ~30
      max.edge = c(5, 20),
      max.n = c(400, 100),
      min.angle = c(30, 21),
      cutoff = 5
    )
  } else if (resolution == "medium") {
    mesh <- INLA::inla.mesh.2d(loc,
      boundary = boundary,
      offset = c(0.0, 30.0),
      ## Shortest correlation range should be ~30
      max.edge = c(5, 20),
      max.n = c(1000, 250),
      min.angle = c(30, 21),
      cutoff = 2
    )
  } else if (resolution == "fine") {
    mesh <- INLA::inla.mesh.2d(loc,
      boundary = boundary,
      offset = c(0.0, 30.0),
      ## Shortest correlation range should be ~30
      max.edge = c(5, 20),
      max.n = c(4000, 500),
      min.angle = c(30, 21),
      cutoff = 1
    )
  } else {
    error("resolution must be \"coarse\", \"medium\", or \"fine\"")
  }
  mesh
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

##' Generate FEM matrices for anisotropic SPDE. Code adapted from TMB anisotropy
##' example.
##'
##' @title Generate anisotropic FEM matrices
##' @param mesh INLA mesh
##' @return List of data expected for \code{spde_aniso_t} TMB object
##' @author John Best
##' @export
generate_aniso_fem <- function(mesh) {
  inla_spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

  ## whi
  dset <- 1:2

  # Triangle info
  TV <- mesh$graph$tv # Triangle to vertex indexing
  V0 <- mesh$loc[TV[, 1], 1:2] # V = vertices for each triangle
  V1 <- mesh$loc[TV[, 2], 1:2]
  V2 <- mesh$loc[TV[, 3], 1:2]
  E0 <- V2 - V1 # E = edge for each triangle
  E1 <- V0 - V2
  E2 <- V1 - V0

  # Calculate Areas
  TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
  Tri_Area <- rep(NA, nrow(E0))
  for (i in 1:length(Tri_Area)) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2 # T = area of each triangle

  G0_inv_diag <- vapply(
    seq_len(nrow(inla_spde$param.inla$M0)),
    function(i) {
      1 / inla_spde$param.inla$M0[i, i]
    }, 0.0
  )
  G0_inv <- as(diag(G0_inv_diag), "TsparseMatrix")

  list(
    n_s      = inla_spde$n.spde,
    n_tri    = nrow(TV),
    Tri_Area = Tri_Area,
    E0       = E0,
    E1       = E1,
    E2       = E2,
    TV       = TV - 1,
    G0       = inla_spde$param.inla$M0,
    ## G0_inv   = as(diag(1 / diag(inla_spde$param.inla$M0)), "TsparseMatrix")
    G0_inv   = G0_inv
  )
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
  if (zero) {
    return(generate_empty_projection(mesh, data_df, group))
  }
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
