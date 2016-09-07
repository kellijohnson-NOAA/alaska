#' Simulate data from a Gompertz population dynamics model
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
#' @param slope
#' @param weightvals
#'
Sim_Gompertz_Fn <- function(n_years, n_stations = 100, phi = NULL,
  SpatialScale = 0.1, SD_O = 0.5, SD_E = 1.0, SD_extra = 1.0,
  rho = 0.5, logMeanDens = 1, Loc = NULL, projection = akCRS,
  slope = 0.5, weightvals = c(10, 2)) {

###############################################################################
## Parameters
###############################################################################
  # Save the input list to access later
  input <- list("phi" = phi, "SpatialScale" = SpatialScale,
    "SD_O" = SD_O, "SD_E" = SD_E, "SD_extra" = SD_extra, "rho" = rho,
    "logMeanDens" = logMeanDens, "projection" = projection,
    "slope" = slope, "weightvals" = weightvals)

  # Determine the starting position from equilibrium
  if (is.null(phi)) phi <- rnorm(1, mean = 0, sd = 1)

  # Calculate the mean growth rate for each subpopulation
  # If all values are the same, then condense alpha to the first value
  alpha <- logMeanDens * (1 - rho)
  if (length(unique(alpha)) == 1) alpha <- alpha[1]

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

  # scale determines the distance at which correlation declines to ~10% of
  # the maximum observed correlation
  # Estimates of "Range" should scale linearly with scale because
  # Range = sqrt(8)/exp(logkappa)
  model_O <- RandomFields::RMgauss(var = SD_O^2, scale = SpatialScale)
  model_E <- RandomFields::RMgauss(var = SD_E^2, scale = SpatialScale)

  # Simulate Omega to obtain an estimate of spatial variation for each location
  # todo: May need to only supply locations that are in a single alpha region
  #       such that Omega is not correlated across boundaries
  RandomFields::RFoptions(spConform = FALSE)
  Omega <- RandomFields::RFsimulate(model = model_O,
    x = Loc[, "x"], y = Loc[, "y"])

  # Simulate Epsilon
  Epsilon <- array(NA, dim = c(n_stations, n_years))
  for(t in 1:n_years) {
    Epsilon[, t] <- RandomFields::RFsimulate(
      model = model_E,
      x = Loc[, "x"], y = Loc[, "y"])
  }
  RandomFields::RFoptions(spConform = TRUE)

  # Determine which subpopulation each location belongs to
  # Find the outer boundaries
  cuts <- unlist(attributes(extent(pol_studyarea))[c("xmin", "xmax")])
  latlimits <- unlist(attributes(extent(pol_studyarea))[c("ymin", "ymax")]) * c(0.0, 1.8)
  # Use quantile to cut into one region per alpha value
  cuts <- floor(quantile(cuts, seq(0, 1, length.out = length(alpha) + 1)))
  # Remove the first value because it represents the lower Longitude limit
  # independent of the number of alpha values
  cuts <- cuts[-1]
  # Create spatial lines for each cut
  if (length(cuts) > 1) {
    # Remove the last value because it represents the higher Longitude limit
    cuts <- cuts[-length(cuts)]
    lines <- sp::SpatialLines(lapply(cuts, function(x) {
      angle <- sample(c(1, runif(1, min = 1 - slope, max = 1 + slope)), 2)
      Lines(Line(cbind(x * angle, latlimits)),
        ID = parent.frame()$i[])
      }))
    proj4string(lines) <- projection
    # create a very thin polygon for each intersected line
    blpi <- rgeos::gBuffer(gIntersection(pol_studyarea, lines, byid = TRUE),
      byid = TRUE, width = 0.000001)
    proj4string(blpi) <- projection
    # split pol_studyarea with each thin polygon
    dpi <- rgeos::gDifference(pol_studyarea, blpi,
      byid = FALSE, drop_lower_td = TRUE)
    proj4string(dpi) <- projection
    # Determine which polygon each point is in
    group <- over(points, sp::disaggregate(dpi))
    pol_grouptrue <- sp::disaggregate(dpi)
  } else {
      group <- rep(1, length.out = NROW(Loc))
      pol_grouptrue <- NULL
  }

###############################################################################
## Calculate Psi
###############################################################################
  Theta <- array(NA, dim = c(n_stations, n_years))
  for (s in 1:n_stations) {
  for(t in 1:n_years) {
    if(t == 1) Theta[s, t] <- as.numeric(
      phi + Epsilon[s, t] + (alpha[group[s]] + Omega[s])/(1 - rho)
      )
    if(t >= 2) Theta[s, t] <- as.numeric(
      rho * Theta[s, t - 1] + alpha[group[s]] + Omega[s] + Epsilon[s, t]
      )
  }}

  # Simulate data
  DF <- array(NA, dim = c(n_stations * n_years, 3),
    dimnames = list(NULL, c("Site", "Year", "Simulated_example")))
  for(s in 1:n_stations) {
  for(t in 1:n_years) {
    counter <- ifelse(s == 1 & t == 1, 1, counter + 1)
    DF[counter, ] <- c(s, t,
      rpois(1, lambda = exp(Theta[s, t] + SD_extra*rnorm(1))))
  }}

  DF <- as.data.frame(DF)
  DF$Simulated_weight <- sapply(DF[, "Simulated_example"], function(x) {
    temp <- rlnorm(x, mean = weightvals[1], sd = weightvals[2])
    temp <- log(temp * exp(weightvals[2]^2 / 2))
    return(sum(temp))
  })
  DF$Longitude <- Loc[DF[, "Site"], 1]
  DF$Latitude <- Loc[DF[, "Site"], 2]

  # Return stuff
  Sim_List <- list("DF" = DF, "phi" = phi, "Loc" = Loc,
    "Omega" = Omega, "Epsilon" = Epsilon, "Theta" = Theta,
    "alpha" = alpha, "cuts" = cuts, "group" = group,
    "pol_studyarea" = pol_studyarea, "pol_grouptrue" = pol_grouptrue,
    "input" = input)

  return(Sim_List)
}
