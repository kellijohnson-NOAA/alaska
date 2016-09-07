#' Create a new polygon based on an existing polygon and a percent increase
#'
#' @param spgeom A polygon
#' @param ratio The change in size, scaled to one, where greater than one is an
#' increase and less than one is a decrease.
#' @param precision The degree you would like your results to be on target
#' @param maxsteps The maximum number of iterations you would like the algorithm
#' to perform
calc_areabuffer <- function(spgeom, ratio = 1.1, precision = 0.0001,
  maxsteps = 20) {
    spgeom.area <- rgeos::gArea(spgeom)
    spgeom.perimeter <- rgeos::gLength(spgeom)
    buffer.width <- (ratio - 1) * spgeom.area / spgeom.perimeter
    buffered.spgeom <- rgeos::gBuffer(spgeom, width = buffer.width)
    achieved.precision <- rgeos::gArea(buffered.spgeom) / spgeom.area - ratio
    steps <- 1
    while(abs(achieved.precision) > precision & steps < maxsteps) {
      buffered.perimeter <- rgeos::gLength(buffered.spgeom)
      buffer.width <- buffer.width - achieved.precision *
        spgeom.area / buffered.perimeter
      buffered.spgeom <- rgeos::gBuffer(spgeom, width = buffer.width)
      achieved.precision <- rgeos::gArea(buffered.spgeom) / spgeom.area - ratio
      steps <- steps + 1
    }
    if(steps == maxsteps & abs(achieved.precision) > precision) {
      warning("Required precision not reached after ", maxsteps,
        " steps. Achieved precision = ", achieved.precision)
    }
    return(buffered.spgeom)
}
