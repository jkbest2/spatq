##' Extract the covariance matrix from an \code{\link[TMB]{sdreport}} and
##' convert it to a correlation matrix.
##'
##' @title Correlation matrix of the fixed-effect parameters
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @return The correlation matrix of the fixed-effect parameters.
##' @author John K Best
##' @importFrom stats cov2cor
##' @export
fixpar_corr <- function(sdr) {
  cov2cor(sdr$cov.fixed)
}

##' Return a data frame with parameter pairs that have a higher absolute
##' correlation than a given cutoff.
##'
##' Note that the returned data frame will not include duplicates. Parameters
##' with repeated names are disambiguated using \code{\link[base]{make.unique}},
##' so they are essentially zero-indexed.
##'
##' @title Fixed-effect parameters with high correlations
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @param corr_cutoff Return pairs of parameters with correlations larger than
##'   this cutoff
##' @return A data frame with columns \code{par1}, \code{par2}, and \code{cor}
##' @author John K Best
##' @export
fixpar_highcorr <- function(sdr, corr_cutoff = 0.8) {
  if (corr_cutoff < 0 || corr_cutoff > 1)
    stop("Cutoff must be between 0 and 1.")
  corr <- fixpar_corr(sdr)
  corr[lower.tri(corr, diag = TRUE)] <- 0
  corr_high <- abs(corr) > corr_cutoff
  indices <- which(corr_high, arr.ind = TRUE)
  parnames <- make.unique(names(sdr$par.fixed))
  df <- apply(indices, 1,
              function(i)
                data.frame(par1 = parnames[i[2]],
                           par2 = parnames[i[1]],
                           corr = corr[i[1], i[2]]))
  df <- Reduce(rbind, df)
  return(df)
}

##' Larger values indicate that this matrix is "more singular", and so more
##' susceptible to numerical error when inverted.
##'
##' @title Calculate the condition number of the fixed-effect parameter
##'   covariance matrix
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @param ... Additional arguments to pass to \code{\link[base]{kappa}}
##' @return The condition number of fixed-effect parameter covariance matrix.
##' @author John K Best
##' @export
fixpar_covcond <- function(sdr, ...) {
  kappa(sdr$cov.fixed, ...)
}

##' Return the square root of the diagonal elements of the fixed-effect
##' parameter covariance matrix.
##'
##' @title Return the standard deviations of the fixed-effect parameters
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @return The named vector of standard deviations
##' @author John K Best
## @export
fixpar_sd <- function(sdr) {
  sds <- sqrt(diag(sdr$cov.fixed))
  names(sds) <- names(sdr$par.fixed)
  sds
}
##' Get the gradients at the estimated values of the fixed effects.
##'
##' @title Gradients at the fitted values
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @return Named vector of gradients
##' @author John K Best
##' @export
fixpar_grad <- function(sdr) {
  sdr$gradient.fixed
}

##' Find the largest gradient at the fitted values.
##'
##' @title Maximum gradient component at fitted parameter values
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @param abs Return the maximum absolute gradient or maximum gradient?
##' @return Largest gradient with name of parameter
##' @author John K Best
##' @importFrom grDevices hcl.colors
##' @importFrom graphics axis image layout lines
##' @export
fixpar_mgc <- function(sdr, abs = TRUE) {
  grad <- fixpar_grad(sdr)
  if (abs) grad <- abs(grad)
  idx <- which.max(grad)
  gr <- grad[idx]
  names(gr) <- names(sdr$par.fixed[idx])
  gr
}

##' Plot the fixed-effect parameter correlation matrix, with pretty labels and
##' clear divisions between parameters. High off-diagonal correlations may
##' indicate a difficult-to-fit model. Uses the "Blue-Red 3" color palette.
##'
##' @title Plot the fixed-effect parameter correlation matrix
##' @param sdr An \code{\link[TMB]{sdreport}} from a TMB model fit
##' @return Nothing; produces plot as a side effect
##' @author John K Best
##' @export
fixpar_corrplot <- function(sdr) {
  npar <- length(sdr$par.fixed)
  parnames <- names(sdr$par.fixed)
  parnamerle <- rle(parnames)
  rlecumsum <- cumsum(c(0, parnamerle$lengths))

  cor <- fixpar_corr(sdr)

  layout(matrix(c(1, 2), ncol = 2), widths = c(0.8, 0.2))
  col <- hcl.colors(20, "Blue-Red 3", rev = TRUE)
  image(seq_len(npar), seq_len(npar), cor, zlim = c(-1, 1),
        axes = FALSE, xlab = NA, ylab = NA, asp = 1, col = col,
        mar = c(1, 1, 0, 0), oma = c(0, 0, 0, 0))
  for (idx in seq_along(parnamerle$values)) {
    endpoints <- c(0.5, parnamerle$lengths[idx] + 0.5) + rlecumsum[idx]
    lab <- par_expr(parnamerle$values[idx])
    axis(1, endpoints, tick = TRUE, labels = FALSE, pos = 0.5)
    axis(1, mean(endpoints), tick = FALSE, labels = lab, pos = 0.5, las = 1)
    axis(2, endpoints, tick = TRUE, labels = FALSE, pos = 0.5)
    axis(2, mean(endpoints), tick = FALSE, labels = lab, pos = 0.5, las = 1)
  }
  lims <- c(0.5, npar + 0.5)
  for (s in (rlecumsum + 0.5)) {
    lines(lims, c(s, s))
    lines(c(s, s), lims)
  }
  plot_corcolorbar(col)
}

##' Plot a colorbar given limits and a set of colors.
##'
##' @title Plot a colorbar legend
##' @param zlim Limits on the color scale
##' @param col Vector of colors
##' @param iasp Inverse aspect ratio (larger values will result in wider colorbars)
##' @author John K Best
##' @export
plot_colorbar <- function(zlim, col = hcl.colors(12), iasp = 1) {
  n <- length(col)
  ## Fix the width and height of each cell to 1; label manually
  x <- c(-0.5, 0.5)
  y <- seq_len(n + 1)
  z <- matrix(seq(min(zlim), max(zlim), length.out = n), nrow = 1)

  zlabs <- signif(seq(min(zlim), max(zlim), length.out = n + 1), 2)

  image(x = x, y = y, z = z,
        col = col, asp = 1 / iasp,
        axes = FALSE, xlab = NA, ylab = NA,
        mar = c(0, 0, 0, 1), oma = c(0, 0, 0, 0))
  axis(4, at = y, labels = zlabs,
       pos = 0.5, lwd = 0, lwd.ticks = 1, las = 1)
}

##' @describeIn plot_colorbar Plot a colorbar for correlations
##' @export
plot_corcolorbar <- function(col = hcl.colors(20, "Blue-Red 3", rev = TRUE)) {
  plot_colorbar(c(-1, 1), col)
}

##' Take the eigendecomposition of the fixed-effect covariance matrix, filter
##' out the smallest eigenvalues, and then order the parameter names by the
##' loadings in the corresponding eigenvectors.
##'
##' The \code{reltol} argument sets a minimum for the ratio of each eigenvalue
##' to the largest eigenvalue. The \code{abstol} argument sets an absolute
##' minimum value. The latter takes precedence if it is passed.
##' @title Find the parameters that contribute to a singular covariance matrix
##' @param x A matrix
##' @param sdr An \code{\link[TMB]{sdreport}} object from a TMB model
##' @param reltol Relative tolerance for small eigenvalues
##' @param abstol Absolute tolerance for small eigenvalues; takes precedence
##' @param vectrans Transformation function applied before sorting eigenvector
##'   elements, usually either `abs` or `identity`
##' @return A list with \code{values} (eigenvalues), \code{vectors}
##'   (eigenvectors), and \code{sortedpars}, which is a matrix of parameter
##'   names with each column ordered by the corresponding loadings of the
##'   eigenvectors (largest loadings at the top)
##' @author John K Best
##' @export
nulleigs <- function(x, reltol = 1e-4, abstol = NULL, vectrans = abs) {
   eigs <- eigen(x)
   if (is.numeric(abstol)) {
     val_idx <- which(eigs$values < abstol)
   } else {
     val_idx <- which((eigs$values / max(eigs$values)) < reltol)
   }
   if (length(val_idx) == 0) {
     message("All eigenvalues above tolerance")
     return(list())
   }

   vecs <- eigs$vectors[, val_idx, drop = FALSE]
   row.names(vecs) <- row.names(x)
   transvecs <- vectrans(vecs)

   ## Use negative to order largest elements first
   ordperms <- apply(-transvecs, 2, order)
   parnames <- make.unique(row.names(x))
   pars <- apply(ordperms, 2, function(r) parnames[r])

   list(values = eigs$values[val_idx], vectors = eigs$vectors[, val_idx],
        sortedpars = pars)
}

##' @describeIn nulleigs Find (nearly) null eigenvalues in the fixed parameter
##'   covariance matrix
##' @export
fixpar_nulleigs <- function(sdr, reltol = 1e-4, abstol = NULL) {
  fpcov <- sdr$cov.fixed
  if (is.null(fpcov)) {
    stop("No `cov.fixed` in `sdr`")
  }
  nulleigs(fpcov, reltol, abstol)
}

##' @describeIn nulleigs Find (nearly) null eigenvalues in the joint precision
##'   matrix
##' @export
jointprec_nulleigs <- function(sdr, reltol = 1e-4, abstol = NULL) {
  jp <- sdr$jointPrecision
  if (is.null(jp)) {
    stop("Error: no joint precision matrix.
          Re-run `sdreport` with `jointPrecision = TRUE`.")
  }
  nulleigs(sdr$jointPrecision, reltol, abstol)
}
