#' @param object An \code{inla.mesh} object
#'
#' @param projection A spatial projection

findlocal <- function(object, projection = akCRS) {

  if (class(object) != "inla.mesh")  stop("object must be of",
    "the class inla.mesh")

  points <- data.frame("x" = object$loc[, 1],
    "y" = object$loc[, 2])
  coordinates(points) <- ~ x + y
  proj4string(points) <- projection

  # Find points in the main part of the mesh
  if (NROW(object$segm$int$idx) > 0) {
    main <- points[object$segm$int$idx[, 1], ]
  } else {
    main <- points[-object$segm$bnd$idx[, 1], ]
  }

  # Create a polygon of the main points
  poly <- SpatialPolygons(list(
    Polygons(list(Polygon(list(main@coords), hole = FALSE)),
    "bound")))
  proj4string(poly) <- proj4string(main)

  # Find which points are inside polygon
  return(ifelse(is.na(over(points, poly)), FALSE, TRUE))
}
