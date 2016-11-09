#' Simulate data from a Gompertz population dynamics model.
#'
#' @description Simulate data for \code{n_years} and \code{n_stations}
#' using the \code{\pkg{RandomFields}} package and a Gompertz population
#' dynamics model. The model uses a recursive equation to simulate population
#' dynamics rather than an autoregressive model. Please ensure that your
#' estimation model uses the same code, otherwise you will estimate biased
#' variance parameters.
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
#' @param rho Density-dependence
#' @param logMeanDens A scalar or vector of log mean density that will be
#' converted into a scalar or vector value for \code{alpha} or mean
#' productivity. \code{alpha} = \code{logMeanDens} * (1 - \code{rho}).
#' @param Loc A two-column matrix of locations.
#' @param projection The projection for your \code{Loc}.
#'
Sim_Gompertz_Fn <- function(n_years, n_stations = 100, phi = NULL,
  SpatialScale = 0.1, SD_O = 0.5, SD_E = 1.0, SD_obs = 1.0,
  rho = 0.5, logMeanDens = 1, Loc = NULL, projection = NULL,
  seed = 1) {

###############################################################################
## Parameters
###############################################################################
  set.seed(seed)
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
  if (!is.null(projection)) sp::proj4string(pol_studyarea) <- projection
  # Find a new polygon that is 10% bigger
  pol_studyarea <- calc_areabuffer(pol_studyarea, ratio = 1.1,
    precision = 0.001)

  # Create SpatialPoints from Loc data
  points <- as.data.frame(Loc)
  coordinates(points) <- ~ x + y
  if (!is.null(projection)) sp::proj4string(points) <- projection

  # Determine which subpopulation each location belongs to
  # Find the outer boundaries
  lonlimits <- unlist(attributes(raster::extent(pol_studyarea))[c("xmin", "xmax")])
  # Make them a little smaller to decrease the likelihood that it will be
  # split bad
  lonlimits[1] <- ifelse(lonlimits[1] < 0 ,
    lonlimits[1] * 0.95, lonlimits[1] * 1.05)
  lonlimits[2] <- ifelse(lonlimits[2] >= 0 ,
    lonlimits[2] * 0.95, lonlimits[2] * 1.05)
  latlimits <- unlist(attributes(raster::extent(
    calc_areabuffer(pol_studyarea, ratio = 3.5)))[c("ymin", "ymax")])
  cuts <- NULL
  if (length(alpha) > 1) {
    table <- 0
    while (any(table < 0.25) | all(table == 1)) {
      cuts <- runif(length(alpha) - 1, min = lonlimits[1], max = lonlimits[2])
      lines_grouptrue <- sp::SpatialLines(lapply(cuts, function(x) {
        Lines(Line(cbind(x, latlimits)),
          ID = parent.frame()$i[])
        }))
      if (!is.null(projection)) sp::proj4string(lines_grouptrue) <- projection
      # Determine which polygon each point is in
      group <- sp::over(points, calc_polys(pol_studyarea, lines_grouptrue))
      table <- table(group) / length(group)
    }
  } else {
      group <- rep(1, length.out = NROW(Loc))
      lines_grouptrue <- NULL
  }
  cuts <- c(latlimits[1], cuts)
  # scale determines the distance at which correlation declines to ~10% of
  # the maximum observed correlation
  # Estimates of "Range" should scale linearly with scale because
  # Range = sqrt(8)/exp(logkappa)
  model_O <- RandomFields::RMgauss(var = SD_O^2, scale = SpatialScale)
  model_E <- RandomFields::RMgauss(var = SD_E^2, scale = SpatialScale)

  RandomFields::RFoptions(spConform = FALSE)
  # Simulate Epsilon
  Epsilon <- array(NA, dim = c(n_stations, n_years))
  for(t in 1:n_years) {
    Epsilon[, t] <- RandomFields::RFsimulate(
      model = model_E,
      x = Loc[, "x"], y = Loc[, "y"])
  }
  Omega <- unlist(split(unlist(tapply(1:NROW(Loc), group,
    function(x) {
    RandomFields::RFsimulate(model = model_O,
      x = Loc[x, "x"], y = Loc[x, "y"])
  })), group))
  RandomFields::RFoptions(spConform = TRUE)

###############################################################################
## Calculate Psi
###############################################################################
  Theta <- array(NA, dim = c(n_stations, n_years))
  DF <- array(NA, dim = c(n_stations * n_years, 9),
    dimnames = list(NULL, c(
      "Site",
      "Year",
      "lambda",
      "group",
      "Epsilon",
      "Omega",
      "alpha",
      "Longitude",
      "Latitude")))
  for (it_s in 1:n_stations) {
  for (t in 1:n_years) {
    if(t == 1) Theta[it_s, t] <- as.numeric(
      phi + Epsilon[it_s, t] + (alpha[group[it_s]] + Omega[it_s])/(1 - rho)
      )
    if(t >= 2) Theta[it_s, t] <- as.numeric(
      rho * Theta[it_s, t - 1] + alpha[group[it_s]] + Omega[it_s] + Epsilon[it_s, t]
      )
    counter <- ifelse(it_s == 1 & t == 1, 1, counter + 1)
    DF[counter, "Site"] <- it_s
    DF[counter, "Year"] <- t
    DF[counter, "lambda"] <- exp(Theta[it_s, t])
    DF[counter, "group"] <- as.numeric(group[it_s])
    DF[counter, "Epsilon"] <- Epsilon[it_s, t]
    DF[counter, "Omega"] <- Omega[it_s]
    DF[counter, "alpha"] <- as.numeric(alpha[group[it_s]])
    DF[counter, "Longitude"] <- Loc[it_s, 1]
    DF[counter, "Latitude"] <- Loc[it_s, 2]
  }}

  DF <- as.data.frame(DF)
  DF <- DF[order(DF$group, DF$Site, DF$Year), ]
  DF$Simulated_example <- rpois(NROW(DF), lambda = DF$lambda)
  DF$encounterprob <- 1 - exp(-DF$lambda)
  DF$zeroinflatedlnorm <- ifelse(DF$Simulated_example > 0, 1, 0) *
    rlnorm(NROW(DF),
      meanlog = log(DF$lambda / DF$encounterprob),
      sdlog = SD_obs)
  DF$phi <- phi
  DF$cuts <- cuts[DF$group]
  DF$sd_O <- SD_O
  DF$sd_E <- SD_E
  DF$sd_obs <- SD_obs
  DF$SpatialScale <- SpatialScale
  DF$seed <- seed

  # Return stuff
  Sim_List <- list("DF" = DF, "lines_grouptrue" = lines_grouptrue)

  return(DF)
}
