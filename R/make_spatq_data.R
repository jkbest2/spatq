make_spatq_data <- function(catch_obs,
                            catch_locs,
                            X_n = NULL,
                            X_w = NULL,
                            Z_n = NULL,
                            Z_w = NULL,
                            R_n = NULL,
                            R_w = NULL,
                            V_n = NULL,
                            V_w = NULL,
                            mesh) {
  spatq_data <- structure(list(), class = "spatqdata")

  return(list(catch_obs = catch_obs,
              X_n = abundance$X_n,
              X_w = abundance$X_w,
              Z_n = abundance$Z_n
              R_n = catchability$R_n,
              R_w = catchability$R_w,
              ))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title FEM matrices from an INLA mesh
##' @param mesh The INLA mesh to use
##' @return A list of sparse FEM matrices from INLA, renamed for compatibility
##'   with TMB's `spde` object
make_fem <- function(mesh) {
  stopifnot(inherits(mesh, "inla.mesh"))
  ## Use INLA::inla.mesh.fem instead of INLA::inla.spde2.matern as recommended by
  ## Finn Lindgren
  fem <- INLA::inla.mesh.fem(mesh)
  ## Rename elements of `fem` so that they are recognized as a `spde` object in
  ## TMB
  names(fem) <- c("M0", "c1", "M1", "M2", "va", "ta")
  return(fem)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Create a sparse matrix that projects from mesh vertices to
##'   observation locations.
##' @param loc Locations to project to
##' @param mesh INLA mesh
##' @param lgl_idx Optional, logical vector of rows to drop
##' @return A sparse (Matrix::dgCMatrix) projection matrix
make_projection <- function(loc, mesh, lgl_idx = NULL) {
  ## Create the full projection matrix
  A <- INLA::inla.spde.make.A(mesh, loc)
  ## Drop rows if necessary (e.g. for subsetting fishery-dependent observations)
  if (!is.null(lgl_idx)) {
    A <- Matrix::drop0(Matrix::Diagonal(x = lbl_idx) %*% A)
  }
  return(A)
}

make_fixeffs <- function(f, data) {
  model.matrix(f, data)
}

make_raneffs <- function(f, data, form = "iid") {
  mm <- model.matrix(f, data)
  attributes(mm) <- form
  if (form = "iid") {
    attributes(mm, n_pars) <- 1
  } else {
    stop("Only form = \"iid\" implemented")
  }
  mm
}
