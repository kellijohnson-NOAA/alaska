#' Exponential population growth with process and observation error
#'
#' Simulate \code{n.pops} with exponential population growth.
#' Populations are observed with observation error and growth
#' with process error and a pre-determined growth rate unique
#' to each population.
#' @author Kelli Faye Johnson
#' @param n.obs Number of observations per population.
#' @param n.pops Number of populations.
#' @param years Length of time series.
#' @param seed Set the seed prior to the simulation.
#' @param growth Growth rates for each population.
#' A vector of length \code{n.pops}.
#' @param sd.obs Standard deviation of observation error.
#' @param sd.process Standard deviation of process error.
#' @export

simulate_data <- function(n.obs = 10, n.pops = 3, 
                          years = 1:10, seed = 199,
                          growth = c(0.5, -0.5, 0.25),
                          sd.obs = 0.5, sd.process = 0.25) {
    set.seed(seed)
    n.years <- length(years)
    n.data <- n.years * n.pops * n.obs
    data <- as.data.frame(matrix(NA, ncol = 4, nrow = n.data))
    colnames(data) <- c("Y", "year", "pop", "rep")
    X <- matrix(NA, ncol = n.years, nrow = n.pops)
    colnames(X) <- years
    X[, 1] <- rnorm(n.pops)
    for(year in 2:length(years)) {
        X[, year] <- X[, year - 1] + growth + rnorm(n.pops, 0, sd.process) 
    }
    Y <- matrix(NA, ncol = n.years, nrow = n.pops * n.obs)
    start <- seq(1, (n.obs*n.pops), by = n.obs)
    end <- seq(n.obs, (n.obs*n.pops), by = n.obs)
    for(pop in 1:n.pops) {
      for(year in 1:n.years){
          Y[(start[pop]:end[pop]), year] <- rnorm(n.obs, X[pop, year], sd.obs)
        }
      }
    data$Y <- as.vector(Y)
    data$year <- rep(years, each = (n.obs * n.pops))
    data$rep <- rep(1:n.obs, times = (n.pops * max(years)))
    data$pop <- rep(rep(1:n.pops, times = n.obs), each = length(years))
    data <- data[order(data$pop, data$rep, data$year), ]
    invisible(list(data = data, X = X))
 }


