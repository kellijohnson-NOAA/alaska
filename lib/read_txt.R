#' Read the \code{.txt} output file from SaTScan
#'
#' Read the output from a SaTScan run and import the results into
#' an \code{R} object for analysis.
#'
#' @param dir The directory where the file is located. The default directory
#' is \code{getwd()}.
#' @param file The filename that you would like to read the results for.
#' The file name can be specified with or without an extension, where
#' \code{.txt} is the proper extension.
#'
#' @examples
#' # Only works on local machine
#' # read_txt(dir = "y25", file = "sim_1-4-1.txt")
read_txt <- function(dir = getwd(), file) {
  # Make sure the file exists
  if (!grepl("\\.", file)) file <- paste0(file, ".txt")
  if (strsplit(file, "\\.")[[1]][2] != "txt") stop("The wrong",
    " file extension was specified for ", file)
  filename <- file.path(dir, file)
  if (!file.exists(filename)) stop("The file ", filename,
    " does not exist.")

  #1. Read the lines
  readin <- scan(filename, what = "", blank.lines.skip = TRUE, sep = ",",
    quiet = TRUE)
  readin <- readin[
    (grep("SUMMARY OF DATA", readin) + 1):
    (grep("PARAMETER SETTINGS", readin) - 2)
  ]
  #2. Pull out the relevant information
  summary <- readin[1:(grep("CLUSTERS DETECTED", readin) - 2)]

  clusters <- apply(
    data.frame("start" = grep("Location IDs included", readin),
    "stop" = grep("Coordinates", readin)), 1, function(x, y = readin) {
      y[x[1]:(x[2] - 1)]
    })
  # clean up
  if (is.matrix(clusters)) {
    clusters <- split(clusters, rep(1:NCOL(clusters), each = NROW(clusters)))
  }
  clusters <- lapply(clusters, function(x) {
      gsub("[[:digit:]]+\\.Location IDs included\\.: |[[:space:]]|r", "", x)
  })
  clusters <- lapply(clusters, function(x) as.numeric(x[x != ""]))

  coordinates <- apply(
    data.frame("start" = grep("Coordinates", readin),
    "stop" = grep("P-value", readin)), 1, function(x, y = readin) {
      y[x[1]:x[2]]
    })
  # clean up
  coordinates <- apply(coordinates, 2, function(x) {
    x[1] <- paste(x[1], x[2], sep = ",")
    x <- x[-2]
    return(x)
  })

  # significance <-
  #   readin[(grep("significance level", readin) + 1):length(readin)]

  #3. Export the results to R
  return(list(
    "summary" = summary,
    "clusters" = clusters,
    "coordinates" = coordinates
  ))
}
