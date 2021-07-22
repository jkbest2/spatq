##' Projects and plots fields. Can plot just estimated fields, or include
##' generative fields for comparison purposes.
##'
##' @title Plot spatial and spatiotemporal fields
##' @param est Vector or named list of parameter estimates
##' @param gen Named list of generative parameter values (not required)
##' @param mesh FEM mesh used for spatial/spatiotemporal fields
##' @param proj_step Size of projection pixels
##' @param par_regex Search string for parameter names (e.g. "omega" will
##'   include both \code{omega_n} and \code{omega_w})
##' @param colorbar Plot a colorbar?
##' @author John K Best
##' @importFrom graphics mtext par
##' @export
plot_field <- function(est, gen = NULL, mesh = generate_mesh(),
                       proj_step = 1, par_regex = NULL,
                       colorbar = FALSE) {
  ## Gather into a named list if it isn't already
  if (!is.list(est))
    est <- gather_nvec(est)
  ## select only the desired parameters if desired
  if (!is.null(par_regex))
    est <- est[grep(par_regex, names(est))]
  ## Use for labeling and extracting relevant parts of `gen`
  parnames <- names(est)

  ## If generative values included, interleave them into a single list
  parvecs <- interleave_estgen(est, gen)

  ## Create an index df for the projection matrix
  nyears <- length(parvecs[[1]]) %/% mesh$n
  if (length(parvecs[[1]]) %% mesh$n)
    stop("Parameter vector and mesh don't match lengths")
  index_df <- create_index_df(proj_step, T = nyears)
  s1 <- unique(index_df$s1)
  s2 <- unique(index_df$s2)
  if (nyears > 1) {
    grps <- index_df$year
  } else {
    grps <- purrr::rep_along(index_df$year, 1)
  }
  projmat <- generate_projection(mesh, index_df, group = grps)

  ## Project each field
  fld <- purrr::map(parvecs,
                    project_field,
                    proj = projmat,
                    dim = c(length(s1), length(s2), nyears))

  nfld <- length(fld)
  parlaboffset <- 1 / nfld
  parlabcoord <- seq(1, parlaboffset, by = -parlaboffset) - (parlaboffset / 2)
  parlab <- purrr::map(seq_along(fld),
                       ~ par_expr(names(fld)[.x], hat = isest(fld[[.x]])))

  lmat <- layout_mat(nyears, !is.null(gen), length(parnames), colorbar)
  lw <- layout_widths(nyears, !is.null(gen), length(parnames), colorbar)
  layout(lmat, lw)
  par(oma = c(1, 2, 1, 1), mar = c(1, 1, 1, 1))
  ## Iterate along parameter names, which will be single or double depending
  ## whether the generative values are included
  for (p in parnames) {
    col <- proc_colormap(p)
    fld_idx <- which(names(fld) == p)
    if (colorbar) {
      zlim <- c(min(unlist(fld[fld_idx])),
                max(unlist(fld[fld_idx])))
      for (idx in fld_idx) {
        apply(fld[[idx]], 3, field_image, zlim = zlim, col = col)
        mtext(parlab[[idx]], side = 2, outer = TRUE, at = parlabcoord[idx])
      }
      plot_colorbar(zlim, col)
    } else {
      ## FIXME Don't need a zlim argument if no colorbar, but probably a
      ## cleaner way to do this.
      for (idx in fld_idx) {
        apply(fld[[idx]], 3, field_image, col = col)
        mtext(parlab[[idx]], side = 2, outer = TRUE, at = parlabcoord[idx])
      }
    }
  }
}

##' Default colormaps for numbers density (viridis) and weight per group
##' (inferno). These two are chosen based on the regexes \code{"_n$"} and
##' \code{"_w$"} respectively. Otherwise uses plasma.
##'
##' @title Colormaps for numbers density and weight per group fields
##' @param parname String like "omega_n" or "epsilon_w"
##' @param n Number of colors in the map
##' @return Vector of colors, as from \code{\link[grDevices]{hcl.colors}}
##' @author John K Best
##' @export
proc_colormap <- function(parname, n = 12) {
  if (grepl("_n$", parname)) {
    col <- hcl.colors(n, "viridis")
  } else if (grepl("_w$", parname)) {
    col <- hcl.colors(n, "inferno")
  } else {
    col <- hcl.colors(n, "plasma")
  }
  col
}

##' Plot the \code{\link[graphics]{image}} of a random field, without axes and
##' with aspect ratio one.
##'
##' @title Plot a field
##' @param z Field matrix
##' @param col Colormap
##' @param ... Additional arguments to pass to \code{\link[graphics]{image}}
##' @author John K Best
##' @export
field_image <- function(z, col = hcl.colors(12), ...) {
  image(z, axes = FALSE, asp = 1, col = col, ...)
}

##' Extracts the elements of \code{gen} with the names in \code{est}, then
##' alternates them so that the first element of the returned list is the first
##' from \code{est}, the second is the element of \code{gen} with the same name,
##' etc. Also adds the \code{"isgen"} attribute (logical) indicating whether
##' these are generative or estimated values.
##'
##' Not symmetric; only elements with names in \code{est} will be kept. All
##' elements in \code{est} should be included in \code{gen}.
##' @title Alternate estimated and generative parameter vectors
##' @param est Named list of estimated parameters
##' @param gen named list of generative parameter vectors
##' @param checklen Check corresponding lengths are equal?
##' @return A list with alternating elements of \code{est} and \code{gen}
##' @author John K Best
##' @export
interleave_estgen <- function(est, gen = NULL, checklen = TRUE) {
  ## Use purrr::map instead of lapply to preserve names
  est <- purrr::map(est, `isest<-`)
  ## If only using estimates, nothing else to do: return early
  if (is.null(gen)) return(est)

  ## Subset gen and reorder to match est
  gen <- gen[names(est)]
  if (checklen) {
    if (length(est) != length(gen)) {
      stop("gen must include all elements of est")
    }
    estlen <- vapply(est, length, 1)
    genlen <- vapply(gen, length, 1)
    if (any(estlen != genlen)) {
      stop("Estimated and generative parameter vectors must be the same length")
    }
  }

  gen <- purrr::map(gen, `isgen<-`)
  ## Concatenate then reorder
  c(est, gen)[order(c(seq_along(est), seq_along(gen)))]
}

##' Useful for projecting multiple-year random effects to grids for plotting
##' with \code{image}.
##'
##' @title Project random effects, reshape to array
##' @param parvec Parameter vector (coerced to vector)
##' @param proj Projection matrix
##' @param dim Output array dimensions (typically \code{c(s1, s2, nyears)})
##' @return Array of dims \code{dim}
##' @author John K Best
##' @export
project_field <- function(parvec, proj, dim) {
  fld <- array(proj %*% as.vector(parvec), dim = dim)
  ## Keep "isgen" attr, if present
  isgen(fld) <- isgen(parvec)
  fld
}

##' Set or get the status of a parameter vector; is it generative or estimated?
##'
##' @title Is this parameter generative or estimated?
##' @param x Object
##' @param value Logical; is this a generative/estimated parameter?
##' @return \code{x} with attribute \code{"isgen"} or logical
##' @author John K Best
##' @export
isgen <- function(x) {
  attr(x, "isgen")
}
##' @describeIn isgen  Is this object estimated?
##' @export
isest <- function(x) !isgen(x)
##' @describeIn isgen Set whether this object is generative?
##' @export
"isgen<-" <- function(x, value = TRUE) {
  attr(x, "isgen") <- value
  return(x)
}
##' @describeIn isgen Set attribute \code{"isgen"} to \code{!isest}
##' @export
"isest<-" <- function(x, value = TRUE) {
  isgen(x) <- !value
  return(x)
}

## Layout functions
##' Construct a matrix for \code{\link[graphics]{layout}} with parameters in
##' rows and years in columns.
##'
##' @title Construct layout matrix for field plots
##' @param nyears Number of years to plot
##' @param incl_gen Include generative fields?
##' @param npar Number of parameters to plot
##' @param colorbar Include colorbar?
##' @return A layout matrix
##' @author John K Best
layout_mat <- function(nyears,
                       incl_gen = TRUE,
                       npar = 1,
                       colorbar = FALSE) {
  ofs <- layout_indoffset(nyears, incl_gen, npar, colorbar)
  ll <- lapply(ofs, function (o) layout_1mat(nyears, incl_gen, o, colorbar))
  lmat <- Reduce(rbind, ll)
  return(lmat)
}
##' @param idxoffset Largest of the previously used indices; this layout matrix
##'   will start at \code{idxoffset + 1}.
##' @describeIn layout_mat Layout matrix for a single parameter
layout_1mat <- function(nyears,
                        incl_gen = TRUE,
                        idxoffset = 0,
                        colorbar = FALSE) {
  s <- matrix(idxoffset + (1:(nyears * (2 - !incl_gen))),
              ncol = nyears, byrow = TRUE)
  if (colorbar)
    s <- cbind(s, rep(max(s) + 1, nrow(s)))
  return(s)
}

##' @describeIn layout_mat Vector of index offsets for each parameter in field
##'   plot
layout_indoffset <- function(nyears,
                             incl_gen = TRUE,
                             npar = 1,
                             colorbar = FALSE) {
  ## Number of indices for each parameter (assume incl_gens are all or none)
  stp <- nyears * (2 - !incl_gen) + colorbar
  seq(0, by = stp, length.out = npar)
}
##' @describeIn layout_mat Widths of field plots; thinner for colorbar if
##'   present
layout_widths <- function(nyears,
                          incl_gen = TRUE,
                          npar = 1,
                          colorbar = FALSE) {
  w1 <- rep(1, nyears * (2 - !incl_gen))
  if (colorbar)
    append(w1, nyears / 10)
  rep(w1, npar)
}
