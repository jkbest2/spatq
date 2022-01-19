##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot array of population states (simulated "truth")
##' @param popstate 3-d array with years in the third dimension
##' @param colorbar Normalize color scale across years and include a colorbar?
##' @param plot_cols How many years in each row?
##' @return
##' @author John K Best
##' @export
plot_popstate <- function(popstate, col = hcl.colors(12), colorbar = TRUE, plot_cols = 5) {
  if (is.array(popstate)) {
    zlim <- c(min(popstate), max(popstate))
    nyears <- dim(popstate)[3]
    popstate <- lapply(seq_len(nyears),
                       function(i) popstate[, , i])
  } else if (is.list(popstate)) {
    nyears <- length(popstate)
    zmin <- min(vapply(popstate, min, 0))
    zmax <- max(vapply(popstate, max, 0))
    zlim <- c(zmin, zmax)
  }
  plot_rows <- ceiling(nyears / plot_cols)

  lmat <- layout_pop(nyears, plot_cols, colorbar)
  layout(lmat)
  par(# mfrow = c(plot_rows, plot_cols),
      oma = c(1, 2, 1, 1),
      mar = c(1, 1, 1, 1))

  if (colorbar) {
    lapply(popstate, field_image, zlim = zlim, col = col)
    plot_colorbar(zlim, col)
  }
}

##' @export
layout_pop <- function(nyears, plot_cols = 5, colorbar = TRUE) {
  plot_rows <- ceiling(nyears / plot_cols)
  plot_tot <- plot_rows * plot_cols
  lm <- matrix(1:plot_tot, byrow = TRUE, nrow = plot_rows, ncol = plot_cols)
  lm[lm > nyears] <- nyears + 2
  cb <- rep(nyears + 1, plot_rows)
  cbind(lm, cb)
}
