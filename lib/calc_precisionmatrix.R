#' Calculate the precision matrix
#'
#' @param spde
#' @param kappa
#' @param tau
#'
calc_precisionmatrix <- function(mesh, kappa, tau) {
  spde <- INLA::inla.spde2.matern(mesh, alpha = 2)

  Q <- INLA::inla.spde2.precision(spde,
    theta = log(c(tau, kappa)))

  # Covariance
  S <- solve(Q)

  # Correlation
  dd <- ((0:1000) / 1000)
  SS <- diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
  SStheory <- (dd * kappa) *
    besselK(dd * kappa, 1)
  SStheory[1] <- 1

  return(list("Q" = Q,
    "distance" = dd,
    "S" = S,
    "SS" = SS, "SStheory" = SStheory))
}
