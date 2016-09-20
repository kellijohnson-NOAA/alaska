#' Create polygons from one main polygon split by lines
#'
#' @param polygon
#' @param lines
#'
calc_polys <- function(polygon, lines) {

  projection <- proj4string(polygon)
  projection.lines <- proj4string(lines)
  if (projection != projection.lines) stop("The projections",
    projection, projection.lines, "do not match")

  # create a very thin polygon for each intersected line
  blpi <- rgeos::gBuffer(
    rgeos::gIntersection(polygon, lines, byid = TRUE),
    byid = TRUE, width = 0.000001)
  proj4string(blpi) <- projection
  # split polygon with each thin polygon
  dpi <- rgeos::gDifference(polygon, blpi,
    byid = FALSE, drop_lower_td = TRUE)
  proj4string(dpi) <- projection

  return(disaggregate(dpi))

}
