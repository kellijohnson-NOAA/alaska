#' Calculate a distance matrix from a set of covariates and a variable
#'
#' A \code{matrix} with three columns is converted into a distance
#' matrix using the locations and a weighting factor for the
#' spatial points
#'
#' @param data A three column matrix or data frame. The first two
#' columns are locations (projected, and not latitude / longitude)
#' and the third column is a variable of interest.
#' @param weight The weight will be used to either increase or
#' decrease the influcence of the spatial locations on the distance
#' matrix. A single value should be supplied because the correlation
#' is assumed to be isotropic.
#' Ono *et al.* (\href{2015}{http://www.sciencedirect.com/science/article/pii/S0165783615001708})
#' used two different weights, 0.1 and 1.0.
#'
#' @author Kelli Faye Johnson
#'
calc_dist <- function(data, weight) {
  if (!is.matrix(data)) data <- as.matrix(data)
  if (NCOL(data) != 3) stop("data needs 3 columns",
    " where your data only has ", NCOL(data))

  # 1. calculate distance for the locations
  locs <- dist(data[, 1:2])^2

  # 2. calculate distance for the variable
  var <- dist(data[, 3])^2

  # 3. weight the locations
  locs <- locs * weight

  # 4. combine the distance matrices
  out <- sqrt(locs + var)

  return(out)
}
