#' Determine and mitigate for missing years in wide formatted data
#'
#' Determine which years are missing and add columns for them.
#' @author Kelli Faye Johnson
#' @param data A dataframe with columns for each year. 
#' @export

    check_missing_years <- function(data = Y, file.log = NA){
        year.range <- as.numeric(range(colnames(data)))
        year.levels <- seq(year.range[1], year.range[2])

        year.values <- colnames(data)
        year.factor <- factor(year.values, levels = year.levels)
        missing.year <- year.levels[which(!year.levels %in% year.values)]
        check.me <- length(year.range[1]:year.range[2]) == dim(data)[2]
        if(!check.me) {
          data.missing <- matrix(NA, ncol = length(missing.year),
                                 nrow = dim(data)[1], 
                                 dimnames = (list(rownames(data), missing.year)))
          data.new <- cbind(data, data.missing)
          data.new <- data.new[, order(as.numeric(colnames(data.new)))]

          if (is.na(file.log)) {
              warning(paste("Survey data contain missing years:", 
                      paste(missing.year, collapse = " ")))
            } else {
             sink(file.log, append = TRUE)
             cat("Survey data contain missing years:", 
                   paste(missing.year, collapse = " "), "\n")
             sink()              
            }
          invisible(data.new)
      } else{
        if (is.na(file.log)) {
            cat("All years were present from:", year.range, "\n")
          } else {
            sink(file.log, append = TRUE)
            cat("All years were present from:", year.range, "\n")
            sink()
          }
          invisible(data)
      }
    }
