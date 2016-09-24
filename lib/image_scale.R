#' Plot a coloured scale bar horizontally for a raster image.
#'
#' @details Add an image scale bar to a plot along the horizontal axis,
#' specified at any height (relative from zero to one) along with numeric
#' text indicating the value that each colour represents.
#'
#' @author Kelli Faye Johnson
#'
#' @param breaks A vector of breaks that encompass the range of values
#' plotted.
#' @param where.y A vector of two numeric values specifying
#' the location of the scale bar on the y axis,
#' where the y axis ranges from zero to one.
#' @param size.font The font size you wish the numeric labels to be.
#' @param n.labels The number of labels you want on the scale bar.
#' @param n.sigfigs The number of significant figures for each label, or
#' a single value that will be repeated for each label.
#' @param colors A vector of 2 character values specifying the
#' colors to use for the color ramp. The values will be supplied to
#' the function \code{grDevices::colorRampPalette}.
#'
#' @return None
#'
#' @examples
#' library(raster)
#' r <- raster(system.file("external/test.grd", package="raster"))
#' breaks <- seq(r@data@min, r@data@max, length.out = 100)
#' colours <- grDevices::colorRampPalette(c("grey100", "grey10"))
#' colours <- colours(length(breaks) - 1)
#' image(r, breaks = breaks, col = colours)
#' image_scale(breaks = breaks, where.y = c(0.95, 1.1), n.labels = 6,
#'   n.sigfigs = 1)
#' dev.off()
#' rm(r, breaks, colours)
#'
image_scale <- function(breaks, where.y = c(0.7, 0.8), size.font = 1.25,
  n.labels, n.sigfigs = 3, colors = c("grey100", "grey10")) {

  xpd.old <- par()$xpd
  on.exit(par(xpd = xpd.old))

  # Make a blank plot
  par(new = TRUE, xpd = NA)
  limits <- c(0, 1)
  plot(limits, limits, type = "n", axes = FALSE, ann = FALSE)

  # Make the colors
  grayscale <- grDevices::colorRampPalette(colors)

  # Make the raster image of the scale bar
  rasterImage(as.raster(matrix(
      grayscale(length(breaks) - 1), nrow = 1)),
      par("usr")[1], par("usr")[4] * where.y[1],
      par("usr")[2], par("usr")[4] * where.y[2])

  # Plot the text of values on the scale image
  labels <- breaks[seq(1, length(breaks), length.out = n.labels)]
  if (length(n.sigfigs) == 1) {
    labels <- round(labels, digits = n.sigfigs)
  } else {
    if (length(labels) != length(n.sigfigs)) stop("n.sigfigs must be",
      " of length equal to the value specified \nin n.labels, or a single",
      " digit to be used for each label.")
    for (isigs in 1:length(labels)) {
      labels[isigs] <- round(labels[isigs], n.sigfigs[isigs])
    }
  }

  x.pos <- ifelse(par("usr")[1] > 0, par("usr")[1] * 1.4, par("usr")[1] * 0.6)
  text(
    x = seq(x.pos, par("usr")[2] * 0.95,
      length.out = n.labels),
    y = par("usr")[4] * mean(where.y),
    cex = size.font,
    labels = labels,
    col = rev(grayscale(n.labels)),
    srt = 90)

}
