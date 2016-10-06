#' Create polygons from one main polygon split by lines
#'
#' @param polygon
#' @param lines
#'
calc_polys <- function(polygon, lines) {

  if (!sp::identicalCRS(polygon, lines)) stop("The projections",
    projection(polygon), projection(lines), "do not match")

  # create a very thin polygon for each intersected line
  blpi <- rgeos::gBuffer(
    rgeos::gIntersection(polygon, lines, byid = TRUE),
    byid = TRUE, width = 0.000001)
  projection(blpi) <- CRS(projection(polygon))
  # split polygon with each thin polygon
  dpi <- rgeos::gDifference(polygon, blpi,
    byid = FALSE, drop_lower_td = TRUE)
  projection(dpi) <- sp::CRS(projection(polygon))

  return(disaggregate(dpi))

}
