#' Simulate data from a Gompertz population dynamics model.
#'
#' @description Simulate data for \code{n_years} and \code{n_stations}
#' using the \code{\pkg{RandomFields}} package and a Gompertz population
#' dynamics model.
#'
#' @details The code for this function originally came from James Thorson
#' and his github repository
#' \url{2016_Spatio-temporal_models/Week 7 -- spatiotemporal models/Lab/Sim_Gompertz_Fn.R}
#'
#' @param n_years The number of years you want data for
#' @param n_stations The number of stations you want samples for
#' @param phi The fraction from equilibrium you want to start from
#' The default is to start at a random \code{rnorm(1, mean = 0, sd = 1)}
#' start value.
#' @param SpatialScale The scale of the spatial random effects, must be
#' in the same units as the locations
#' @param SD_O The marginal variance of Omega.
#' @param SD_E The marginal variance of temporal and spatial process error.
#' @param SD_extra
#' @param rho Density-dependence
#' @param logMeanDens
#' @param Loc
#' @param projection
#'
Sim_Gompertz_Fn <- function(n_years, n_stations = 100, phi = NULL,
  SpatialScale = 0.1, SD_O = 0.5, SD_E = 1.0, SD_extra = 1.0, SD_obs = 1.0,
  rho = 0.5, logMeanDens = 1, Loc = NULL, projection = akCRS) {

###############################################################################
## Parameters
###############################################################################
  # Save the input list to access later
  input <- list("phi" = phi, "SpatialScale" = SpatialScale,
    "SD_O" = SD_O, "SD_E" = SD_E, "SD_extra" = SD_extra, "SD_obs" = SD_obs,
    "rho" = rho,
    "logMeanDens" = logMeanDens, "projection" = projection)

  # Determine the starting position from equilibrium
  if (is.null(phi)) phi <- rnorm(1, mean = 0, sd = 1)

  # Calculate the mean growth rate for each subpopulation
  # If all values are the same, then condense alpha to the first value
  alpha <- logMeanDens * (1 - rho)
  if (length(unique(alpha)) == 1) alpha <- alpha[1]
  if (length(alpha) > 1) {
    percentinc <- (alpha[2] - alpha[1]) / alpha[1] * 100
  } else {
    percentinc <- 0
  }

###############################################################################
## Spatial model
###############################################################################
  # Randomly generate the locations if a matrix is not given
  # such that each polygon is approximately square, and if there is
  # only one alpha value then the spatial landscape is 1 x 1
  if (is.null(Loc)) {
    Loc <- cbind(
      "x" = runif(n_stations, min = 0, max = 1),
      "y" = runif(n_stations, min = 0, max = 1 / length(alpha)))
  } else {
    # If locations are given, determine how many stations and
    # set column names
    n_stations <- NROW(Loc)
    if (NCOL(Loc) != 2) stop("Loc does not have two columns")
    colnames(Loc) <- c("x", "y")
  }
  # Create a polygon with a buffer around the locations
  pol_studyarea <- as(raster::extent(Loc), "SpatialPolygons")
  proj4string(pol_studyarea) <- projection
  # Find a new polygon that is 10% bigger
  pol_studyarea <- calc_areabuffer(pol_studyarea, ratio = 1.1,
    precision = 0.001)

  # Create SpatialPoints from Loc data
  points <- as.data.frame(Loc)
  coordinates(points) <- ~ x + y
  proj4string(points) <- projection

  # Determine which subpopulation each location belongs to
  # Find the outer boundaries
  lonlimits <- unlist(attributes(extent(pol_studyarea))[c("xmin", "xmax")])
  latlimits <- unlist(attributes(extent(
    calc_areabuffer(pol_studyarea, ratio = 3.5)))[c("ymin", "ymax")])
  if (length(alpha) > 1) {
    table <- 0
    while (any(table < 0.25)) {
      cuts <- runif(length(alpha) - 1, min = lonlimits[1], max = lonlimits[2])
      lines_grouptrue <- sp::SpatialLines(lapply(cuts, function(x) {
        Lines(Line(cbind(x, latlimits)),
          ID = parent.frame()$i[])
        }))
      projection(lines_grouptrue) <- projection
      # Determine which polygon each point is in
      group <- over(points, calc_polys(pol_studyarea, lines_grouptrue))
      table <- table(group) / length(group)
    }
  } else {
      group <- rep(1, length.out = NROW(Loc))
      lines_grouptrue <- NULL
      cuts <- NA
  }

  # scale determines the distance at which correlation declines to ~10% of
  # the maximum observed correlation
  # Estimates of "Range" should scale linearly with scale because
  # Range = sqrt(8)/exp(logkappa)
  model_O <- RandomFields::RMgauss(var = SD_O^2, scale = SpatialScale)
  model_E <- RandomFields::RMgauss(var = SD_E^2, scale = SpatialScale)

  RandomFields::RFoptions(spConform = FALSE)
  Omega <- unlist(split(unlist(tapply(1:NROW(Loc), group,
    function(x, data = Loc, model = model_O, sd = SD_O) {
    RandomFields::RFsimulate(model = model,
      x = data[x, "x"], y = data[x, "y"]) - sd^2/2
  })), group))

  # Simulate Epsilon
  Epsilon <- array(NA, dim = c(n_stations, n_years))
  for(t in 1:n_years) {
    Epsilon[, t] <- RandomFields::RFsimulate(
      model = model_E,
      x = Loc[, "x"], y = Loc[, "y"])
  }
  RandomFields::RFoptions(spConform = TRUE)

###############################################################################
## Calculate Psi
###############################################################################
  Theta <- array(NA, dim = c(n_stations, n_years))
  DF <- array(NA, dim = c(n_stations * n_years, 3),
    dimnames = list(NULL, c("Site", "Year", "lambda")))
  for (s in 1:n_stations) {
  for (t in 1:n_years) {
    if(t == 1) Theta[s, t] <- as.numeric(
      phi + Epsilon[s, t] + (alpha[group[s]] + Omega[s])/(1 - rho)
      )
    if(t >= 2) Theta[s, t] <- as.numeric(
      rho * Theta[s, t - 1] + alpha[group[s]] + Omega[s] + Epsilon[s, t]
      )
    counter <- ifelse(s == 1 & t == 1, 1, counter + 1)
    DF[counter, "Site"] <- s
    DF[counter, "Year"] <- t
    DF[counter, "lambda"] <- exp(Theta[s, t] + SD_extra*rnorm(1))
  }}

  DF <- as.data.frame(DF)

  DF$Simulated_example <- rpois(NROW(DF), lambda = DF$lambda)
  DF$encounterprob <- 1 - exp(-DF$lambda)
  DF$zeroinflatedlnorm <- ifelse(DF$Simulated_example > 0, 1, 0) *
    rlnorm(NROW(DF),
      meanlog = log(DF$lambda / DF$encounterprob),
      sdlog = SD_obs)
  DF$Longitude <- Loc[DF[, "Site"], 1]
  DF$Latitude <- Loc[DF[, "Site"], 2]

  # Return stuff
  Sim_List <- list("DF" = DF, "phi" = phi, "Loc" = Loc,
    "Omega" = Omega, "Epsilon" = Epsilon, "Theta" = Theta,
    "alpha" = alpha, "cuts" = cuts, "group" = group,
    "input" = input, "lines_grouptrue" = lines_grouptrue,
    "n_grouptrue" = length(alpha),
    "percentinc" = percentinc)

  return(Sim_List)
}
