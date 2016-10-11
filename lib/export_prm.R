#' Export a \code{.prm} file with the appropriate file names inserted.
#'
#' The \code{.prm} file must have each input and output file specified
#' using full filepath names. The function \code{export_prm} opens a
#' template file and changes all of the files to the appropriate file
#' names such that the model can be run.
#'
#' @param file_in A full file path to the \code{.prm} file that you
#' would like to use as a template.
#' @param dir The directory where you would like the files saved.
#' @param file A character value specifying the file name to be used,
#' without an extension. The same file name will be used for all extensions.
#' @param shape A character value of either cirlce or ellipse that dictates
#' how the spatial clusters will be defined.
#' @param size A numeric value between 0 and 100 that dictates the
#' maximum spatial size in population at risk. SaTScan reccommends values
#' less than or equal to 50 with 50 being the default.
#'
#' @examples
#' # This example will only work on my local machine.
#' # export_prm("d:/alaska/data/normal.prm", getwd(), "new")
#'
export_prm <- function(file_in, dir, file,
  shape = c("circle", "ellipse"), size = 50) {
  export <- file.path(dir, paste0(file, ".prm"))
  # Change to forward slashes
  dir <- gsub("/", "\\\\", dir)

  if (size < 0 | size > 100) stop("size must be between [0:100]")
  shape <- match.arg(shape, choices = c("circle", "ellipse"))

  #1. Read in the template
  if (!file.exists(file_in)) stop("The file ", file_in, " does not exist.")
  prm <- readLines(con = file_in)
  cas <- grep(pattern = "\\.cas", x = prm)
  geo <- grep(pattern = "\\.geo", x = prm)
  txt <- grep(pattern = "\\.txt", x = prm)
  shp <- grep(pattern = "SpatialWindowShapeType", x = prm)
  sze <- grep(pattern = "MaxSpatialSizeInPopulationAtRisk", x = prm)
  cor <- grep(pattern = "CoordinatesType", x = prm)

  #2. Change the file names of the .geo, .cas, and .txt files
  prm[cas] <- paste0("CaseFile=",
    file.path(dir, paste0(file, ".cas"), fsep = "\\"))
  prm[geo] <- paste0("CoordinatesFile=",
    file.path(dir, paste0(file, ".geo"), fsep = "\\"))
  prm[txt] <- paste0("ResultsFile=",
    file.path(dir, paste0(file, ".txt"), fsep = "\\"))
  prm[shp] <- paste0("SpatialWindowShapeType=",
    ifelse(shape == "circle", 0, 1))
  prm[sze] <- paste0("MaxSpatialSizeInPopulationAtRisk=", size)
  prm[cor] <- paste0("CoordinatesType=", ifelse(shape == "circle", 1, 0))

  #3. Export the new .prm file with the appropriate name.
  writeLines(prm, export)
}
