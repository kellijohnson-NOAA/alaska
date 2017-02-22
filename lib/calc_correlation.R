#' Calculate covariances and correlations of a sparse matrix
#'
#' To do: more here.
#'
#' @param mesh
#' @param kappa
#' @param tau
#' @param plot A file path for the plot, if you wish to generate one.
#'
calc_correlation <- function(mesh, kappa, tau, plot = NULL,
  model = NULL) {

  # Construct a Q matrix
  # Supply parameters on the log scale
  Q <- INLA::inla.spde2.precision(
    INLA::inla.spde2.matern(mesh, alpha = 2),
    theta = log(c(tau, kappa)))

  # Perform some fancy reordering
  reo=inla.qreordering(Q, reordering="metis")
  # Need to invert the indexing:
  neworder <- reo$reordering
  neworder[neworder] <- 1:length(neworder)
  # Reorder the matrix:
  Q.reordered <- Q[neworder, neworder]

  # Reference point for covariance/correlation comparisons:
  ref.s <- which.min(
    (mesh$loc[, 1] - mean(range(mesh$loc[, 1])))^2 +
    (mesh$loc[, 2] - mean(range(mesh$loc[, 2])))^2)

  # Calculate covariances (S) and correlations (SS):
  S = solve(Q)
  SS = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
  D = as.matrix(dist(mesh$loc))

  # Theoretical Matern correlations and covariances:
  dd <- (0:1000) / 1000
  SS.theory <- (dd * kappa) * besselK(dd * kappa, 1)
  SS.theory[1] <- 1
  S.theory <- SS.theory /
    (4 * pi * kappa^2) / tau^2

  if (!is.null(plot)) {
    if (!grepl("\\.png$", plot)) {
      warning("Adding the .png extension onto ", plot)
      plot <- paste0(plot, ".png")
    }
    png(file = plot)
    plot(D[ref.s, ], SS[ref.s, ], type = "p", pch = 20,
      xlab = "Distance", ylab = "Correlation",
      ylim = c(-0.005, 1.005))
    lines(dd, SS.theory, type = "l",
      col = rgb(0.5, 0.5, 0.5), lwd = 2)
    if (!is.null(model)) {
      lines(model, col = "red")
    }
    dev.off()
  }

  return(list("covariance" = S, correlation = SS))
}
