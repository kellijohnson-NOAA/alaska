#' Split outer limits by a line at \code{cuts}.
#'
#' @param cuts Where you want the lines.
#' @param limits The outer limits, in which to place the lines
#' @param projection The projection of your data.
#'
calc_lines <- function(cuts, limits, projection = NULL) {
  lines <- list()
  for (ii in seq_along(cuts)) {
    lines[[ii]] <- sp::Lines(sp::Line(cbind(cuts[ii], limits)), ID = ii)
  }
  lines <- sp::SpatialLines(lines)
  if (!is.null(projection)) sp::proj4string(lines) <- projection
  return(lines)
}
