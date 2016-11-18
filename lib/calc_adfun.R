#' Make a tmb object to fit using optim or some other algorithm in R.
#'
#' @param data A data set generated from \code{\link{Sim_Gompertz_Fn}}.
#' @param mesh A triangulation network generated from \code{\pkg{INLA}}.
#' @param tmb The library name that was loaded using
#' \code{\link[TMB]{dynlib}}.
#' @param variable A character value of the dependent variable that you
#' want to model. For example, \code{"Simulated_example"}. The character
#' value must be a column name of \code{data}.
#'
calc_adfun <- function(data, mesh, tmb, variable) {
  # Start with data
  ## Make sure that the variable name is in \code{data}
  if(!variable %in% colnames(data)) {
    stop(variable, "is not a column of data")
  }
  ## Make sure that Year is a column name of \code{data}
  if(!"Year" %in% colnames(data)) {
    stop("Year is not a column of data")
  }
  ## Make sure that the mesh is of class inla.mesh
  if(class(mesh) != "inla.mesh") {
    stop("mesh is not of class inla.mesh")
  }
  ## Determine if using counts or weights
  what <- ifelse(all(is.integer(data[, variable])), 0, 1)
  ## Create the spde object from the mesh
  spde <- INLA::inla.spde2.matern(mesh)

  Data <- list(
    Options_vec = what,
    n_i = NROW(data),
    n_x = mesh$n,
    n_t = max(data$Year),
    n_p = 1,
    x_s = mesh$idx$loc - 1,
    c_i = data[, variable],
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
