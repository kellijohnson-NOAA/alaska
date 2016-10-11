#' Export point locations and values to space deliminated files
#'
#' The points must have a unique identifer and be deliminated
#' using spaces or tabs. The points will then be read into
#' SatScan and used for a spatial analysis.
#' The points can be projected or in latitude and longitude,
#' though decimal degrees should be used if they are in the later
#' form. Rownames will be used as identifiers for both the spatial
#' locations and the values of interest.
#'
#' @param points A SpatialPointsDataFrame.
#' @param dir The directory where you would like the files saved.
#' @param file A character value specifying the file name to be used,
#' without an extension. The same file name will be used for both extensions.
#' @param get A character value specifying the name of the variable
#' that should be extracted from the \code{@data} data frame of \code{points}.
#' @param projection The data must be projected using latitude and longitude
#' in decimal degrees for the program to export a shape file.
#'
#' @examples
#' # This example will only work on my local machine.
#' # export_geo(best[[51]]$info, getwd(), "test")
#'
export_geo <- function(points, dir, file, get = "omega",
  projection = sp::CRS("+proj=longlat +ellps=WGS84")) {
  filepath <- file.path(dir, paste0(file, ".geo"))
  points <- sp::spTransform(points, projection)
  data <- points@coords
  val <- 1:2
  if (identical(projection, sp::CRS("+proj=longlat +ellps=WGS84"))) {
    val <- rev(val)
  }
  data <- cbind(paste0("r", 1:NROW(data)), data[, val[1]], data[, val[2]])
  write.table(data, quote = FALSE,
    file = filepath, sep = " ", row.names = FALSE, col.names = FALSE)

  data <- points@data[, get]
  data <- data.frame(paste0("r", seq_along(data)), 1, data)
  filepath <- gsub(".geo", ".cas", filepath)
  write.table(data, quote = FALSE,
    file = filepath, sep = " ", row.names = FALSE, col.names = FALSE)
}
