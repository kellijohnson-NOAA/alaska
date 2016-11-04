#' Make a tmb object to fit using optim or some other algorithm in R.
#'
#' @param data
#' @param mesh
#' @param tmb
#' @param variable A character value of the dependent variable that you
#' want to model. For example, \code{"Simulated_example"}.
#'
calc_adfun <- function(data, mesh, tmb, variable) {
  # Start with data
  ## Determine if using counts or weights
  what <- ifelse(all(is.integer(data[, variable])), 0, 1)
  ## Create the spde object from the mesh
  spde <- inla.spde2.matern(mesh)

  Data <- list(
    Options_vec = what,
    n_i = NROW(data),
    n_x = mesh$n,
    n_t = max(data$Year),
    n_p = 1,
    x_s = mesh$idx$loc - 1,
    c_i = data[, variable],
    s_i = mesh$idx$loc - 1,
    t_i = data[, "Year"] - 1,
    X_xp = matrix(1, ncol = 1, nrow = mesh$n),
    G0 = spde$param.inla$M0,
    G1 = spde$param.inla$M1,
    G2 = spde$param.inla$M2)

  # Calculate the priors
  Parameters <- calc_priors(Data)

  # Determine if using counts or weights
  Map <- NULL
  if (Data$Options_vec == 0) {
    Map[["theta_z"]] <- factor(c(NA, NA))
  }

  # Make object for TMB
  tmbobject <- MakeADFun(data = Data, parameters = Parameters,
    random = c("Epsilon_input", "Omega_input"),
    map = Map, hessian = FALSE, DLL = tmb, silent = TRUE)

  return(tmbobject)
}
