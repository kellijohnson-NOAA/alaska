#' Read in results for alaska simulation
#'
#' @param data Specify a simulated data set, more than likely named
#' \code{sim_data}.
#' @param report Specify the results of running the spatial model,
#' more than likely named \code{Report}.
#' @param dir A directory where the \code{.RData} file specified in
#' \code{file} is stored. If the full path is provided in \code{file}
#' then set \code{dir = NULL}.
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
    if (exists("sim_data", inherits = FALSE)) rm(sim_data)
    if (exists("Report", inherits = FALSE)) rm(Report)

    if (is.null(dir)) { path <- file } else { path <- file.path(dir, file)}
    load(path)
    data <- sim_data
    report <- Report
  }

  # Combine data and report into a single list, such that the function
  # can return everything that is imported.
  # At a minimum the user must supply a \code{report} object.
  if (!is.list(report)) {stop("report must be a list")}
  # Remove the mesh object from data because it is just a replicate of
  # the mesh that is provided in report$mesh
  if (!is.null(data)) data <- data[-which(names(data) == "mesh")]
  results <- append(data, report)
  rm(data); rm(report)

###############################################################################
## Subset the data based on the interior of the mesh
###############################################################################
  # Get the points inside the blue oval within the mesh
  localboundaries <- findlocal(object = results$mesh, projection = projection)

  # Get all information for each mesh location
  info <- data.frame(
    "x" = results$mesh$loc[, 1],
    "y" = results$mesh$loc[, 2],
    "omega" = results[["Omega_x"]],
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
  if (is.null(results$pol_grouptrue)) {
    info@data$true <- 1
  } else {info@data$true <- over(info, results$pol_grouptrue)}

###############################################################################
## Save information to the results list
###############################################################################
  # Create a vector of breaks for the estimated values of Omega
  results$breaks <- seq(min(results[["Omega_x"]]),
    max(results[["Omega_x"]]), length.out = 100)

  results$info <- info
  results$localboundaries <- localboundaries
  names(results)[which(names(results) == "Loc")] <- "loc_true"
  results$file <- file

  return(results)
}
