#' @param object An \code{inla.mesh} object
#'
#' @param projection A spatial projection

findlocal <- function(object, projection = akCRS) {

  if (class(object) != "inla.mesh")  stop("object must be of",
    "the class inla.mesh")

  points <- data.frame("x" = object$loc[, 1],
    "y" = object$loc[, 2])
  coordinates(points) <- ~ x + y
  projection(points) <- projection

  # Create a polygon of the main points
  poly <- calc_meshbound(mesh = object, projection = projection)$poly

  # Find which points are inside polygon
  return(ifelse(is.na(over(points, poly)), FALSE, TRUE))
}
