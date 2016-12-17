#' Read in results for alaska simulation
#'
#' @param report Specify the results of running the spatial model,
#' more than likely named \code{Report}.
#' @param dir A directory where the \code{.RData} file specified in
#' \code{file} is stored. If the full path is provided in \code{file}
#' then set \code{dir = NULL}.
#' @param file The file name of the \code{.RData} object that houses
#' \code{report} for a given simulation replicate.
#' @param projection The spatial projection of the data
#'
read_results <- function(report = NULL,
  dir = getwd(), file = NULL, projection = "NA") {

  results <- list()

###############################################################################
## Determine if the function is loading saved data or using already loaded data
###############################################################################
  if (!is.null(file)) {
    if (exists("Report", inherits = FALSE)) rm(Report)

    if (is.null(dir)) {
      path <- file
    } else {
      path <- file.path(dir, file)
    }
    load(path)
    results <- Report
    rm(Report)
  } else {results <- report; rm(report)}

  # Combine data and report into a single list, such that the function
  # can return everything that is imported.
  # At a minimum the user must supply a \code{report} object.
  if (!is.list(results)) {stop("report must be a list")}

  # Read in the report file, which stores the model output and combine
  # it with the data used to fit the model.
  em <- read_report(results)
  combine <- merge(results$DF, em,
    by = c("Year", "Site"),
    suffixes = c("_om", "_em"),
    all = TRUE)

###############################################################################
## Subset the data based on the interior of the mesh
###############################################################################
  # Get the points inside the blue oval within the mesh
  localboundaries <- data.frame("Site" = 1:results$mesh$n,
    "local" = findlocal(object = results$mesh, projection = projection))
  combine <- merge(combine, localboundaries, by = "Site", all = TRUE)
  info <- data.frame(
    "x" = results$mesh$loc[, 1],
    "y" = results$mesh$loc[, 2],
    "omega" = results[["Omega_x"]])
  # Subset by those inside the blue line
  results$info <- info[localboundaries$local, ]
  sp::coordinates(info) <- ~ x + y
  raster::projection(info) <- projection

###############################################################################
## Save information to the results list
###############################################################################
  combine$n_years <- NCOL(results[["Epsilon_xt"]])

  results$localboundaries <- localboundaries
  names(results)[which(names(results) == "Loc")] <- "loc_true"
  results$file <- file
  results$n_years <- NCOL(results[["Epsilon"]])
  results$all <- combine

return(results)

}
