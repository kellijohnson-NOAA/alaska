#' Read in results for alaska simulation
#'
#' @param data Specify a simulated data set, more than likely named
#' \code{sim_data}.
#' @param report Specify the results of running the spatial model,
#' more than likely named \code{Report}.
#' @param dir A directory where the \code{.RData} file specified in
#' \code{file} is stored.
#' @param file The file name of the \code{.RData} object that houses
#' \code{data} and \code{report} for a given simulation replicate.
#' @param size A numeric value specifying the cell size of the spatial
#' grid. The value will be repeated to create square grid cells.
#' @param dimension A vector of two values, specifying the number of
#' grid cells in the \code{x} and \code{y} dimensions of the spatial grid.
#' @param projection The spatial projection of the data
#'
#' @examples
#' information <- read_results(data = sim_data, report = Report)
#'
read_results <- function(data = NULL, report = NULL,
  dir = getwd(), file = NULL, size = 0.01, dimension = c(130, 50),
  projection = akCRS) {

  results <- list()

###############################################################################
## Determine if the function is loading saved data or using already loaded data
###############################################################################
  if (!is.null(file)) {
    rm(sim_data); rm(Report)
    load(file.path(dir, file))
    data <- sim_data
    report <- Report
  }

###############################################################################
## Subset the data based on the interior of the mesh
###############################################################################
  # Get the points inside the blue oval within the mesh
  localboundaries <- findlocal(object = data$mesh, projection = projection)

  # Get all information for each mesh location
  info <- data.frame(
    "x" = data$mesh$loc[, 1],
    "y" = data$mesh$loc[, 2],
    "omega" = report[["Omega_x"]],
    "clustering" = NA)
  # Subset by those inside the blue line
  info <- info[localboundaries, ]
  coordinates(info) <- ~ x + y
  projection(info) <- projection

  # Create a grid for plotting later
  results$info.grid  <- SpatialGrid(GridTopology(
    cellcentre.offset = bbox(info)[, "min"],
    cellsize = c(size, size), cells.dim = dimension))

  # Calculate the mesh points that are in each true group
  if (is.null(data$pol_grouptrue)) {
    info@data$true <- 1
  } else {info@data$true <- over(info, data$pol_grouptrue)}

###############################################################################
## Save information to the results list
###############################################################################
  # Create a vector of breaks for the estimated values of Omega
  results$breaks <- seq(min(report[["Omega_x"]]),
    max(report[["Omega_x"]]), length.out = 100)

  results$info <- info
  results$alpha <- data$alpha
  results$SigmaO <- report$SigmaO
  results$SigmaE <- report$SigmaE
  results$rho <- report$rho
  results$mesh <- data$mesh
  results$group <- data$group
  results$localboundaries <- localboundaries
  results$loc_true <- data$Loc

  return(results)
}
