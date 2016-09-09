#' Determine the points that creeate the inner boundary for an \code{inla.mesh}
#' object. The boundary is typically delineated by a blue line.
#'
#' @param mesh An \code{inla.mesh} object.
#' @param projection The projection you would like to use for the resulting
#' polygon. The default is \code{NULL}, and if \code{projection = NULL}, then
#' the function will return unprojected points and empty objects for the
#' two types of polygons.
#'
calc_meshbound <- function(mesh, projection = NULL) {

  if (class(mesh) != "inla.mesh") {
    stop("mesh is not an inla.mesh object")
  }

  # Find points in the main part of the mesh
  if (NROW(mesh$segm$int$idx) <= 0) {
    locators <- unique(c(mesh$segm$bnd$idx))
  } else {
    locators <- unique(unlist(mesh$segm$int$idx))
  }

  points <- mesh$loc[locators, ]
  points <- as.data.frame(points)

  colnames(points) <- c("x", "y", "ID")

  if (!is.null(projection)) {
    points$num <- 1:NROW(points)
    coordinates(points) <- ~ x + y
    proj4string(points) <- projection

    polygons <- SpatialPolygons(list(
    Polygons(list(Polygon(list(points@coords))),
    "bound")))

    # Create a convex hull polygon that does not follow the
    # points exactly
    hull <- SpatialPolygons(
      tapply(1:length(points$ID), points$ID,
      function(x, data = points@coords) {
        X <- data[x, ]
        X.chull <- chull(X)
        X.chull <- c(X.chull, X.chull[1])
        Polygons(list(Polygon(X[X.chull, ])), parent.frame()$i[])
    }))
    projection(hull) <- projection
    # Create a polygon that follows the points exactly
    poly <- SpatialPolygons(
      tapply(1:length(points$ID), points$ID,
      function(x, data = points@coords) {
        data <- data[!duplicated(apply(data, 2, paste0)), ]
        Polygons(list(Polygon(data)), parent.frame()$i[])
    }))
    projection(poly) <- projection
  } else {
    hull <- NULL
    poly <- NULL
  }

  return(list("points" = points, "hull" = hull, "poly" = poly))

}
