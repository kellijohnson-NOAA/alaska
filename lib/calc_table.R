#' Calculate an xtable and transform it into a word table.
#'
#' @param table A \code{data.frame}
#' @param file A character value providing the file name with
#' no extension.
#' @param dir A directory where the result should be stored.
#'
#' @examples
#' calc_table(data.frame("yes" = 1:3, "$\\Omega$" = 4:6), "temp", getwd())
#' unlink(dir(pattern = "temp"), recursive = FALSE, force = FALSE)
#'
calc_table <- function(table, file, dir, ...) {
  dir.current <- getwd()
  on.exit(setwd(dir.current))

  if (!file.exists(dir)) stop(dir, " does not exist")
  setwd(dir)

  tex <- paste0(file, ".tex")
  doc <- paste0(file, ".doc")

  print(xtable(table, ...), file = tex,
    include.rownames = FALSE, sanitize.text.function = function(x){x})
  system(paste("pandoc", tex, "-o", doc))
}
