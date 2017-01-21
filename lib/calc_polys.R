#' Create polygons from one main polygon split by lines
#'
#' @param polygon
#' @param lines
#'
calc_polys <- function(polygon, lines) {

  if (!sp::identicalCRS(polygon, lines)) stop("The projections",
    sp::proj4string(polygon), sp::proj4string(lines), "do not match")

  # create a very thin polygon for each intersected line
  blpi <- rgeos::gBuffer(
    rgeos::gIntersection(polygon, lines, byid = TRUE),
    byid = TRUE, width = 0.000001)
  sp::proj4string(blpi) <- CRS(sp::proj4string(polygon))
  # split polygon with each thin polygon
  dpi <- rgeos::gDifference(polygon, blpi,
    byid = FALSE, drop_lower_td = TRUE)
  sp::proj4string(dpi) <- sp::CRS(sp::proj4string(polygon))

  return(disaggregate(dpi))

}
