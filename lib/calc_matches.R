#' Determine the number of matches among two columns.
#'
#' The two columns may have different values such as \code{1:2} and
#' \code{3:4}, but the function will determine which levels align to
#' make the most matches and then assign the estimated names to be the
#' same as the true names. The code only works for numeric values.
#'
#' @param data A data frame that either contains two columns, with the
#' first column being the true value and the second column being the
#' estimated value, or a data frame with many columns, where the columns
#' \code{"true"} and \code{"est"} exist within the data frame.
#'
#' @returns A data frame with all of the original columns and a new column
#' named \code{"estgroup"} providing the new names for the estimates.
#' If no match is found the estimate will be assigned to a new
#' category.
#'
#' @author Kelli Faye Johnson
#'
#' @examples
#' data <- data.frame("true" = rep(1:2, length = 6),
#'   "est" = rep(1:3, each = 4))
#' calc_matches(data)
#'
calc_matches <- function(data) {
  # 01. Check the data has the proper format
  if (is.null(colnames(data))) colnames(data) <- c("true", "est")
  if (!any(c("true", "est") %in% colnames(data))) {
    stop("The data does not contain the appropriate information,\n",
      "the columns 'true' and 'est' must be present.")
  }
  info <- data[, c("true", "est")]
  works <- c("numeric", "integer")
  if (!class(info$true) %in% works | !class(info$est) %in% works) {
    stop("calc_matches only works for numeric values")
  }

  # 02. Find the number of matches
  trues <- unique(info$true); trues <- trues[order(trues)]
  ests <- unique(info$est); ests <- ests[order(ests)]
  match.table <- t(apply(table(info$true, info$est), 1, prop.table))

  # 03. Duplicate the match.table
  matches <- match.table
  colnames(matches) <- rep(0, NCOL(matches))

  # 04. Generate new labels for each match and remove rows from
  # match.table one at a time until all of the highest number of
  # matches are assigned.
  while (!max(match.table, na.rm = TRUE) == 0) {
    max <- which(match.table == max(match.table), arr.ind = TRUE)[1, ]
    colnames(matches)[as.numeric(colnames(match.table)[max[2]])] <-
      rownames(match.table)[as.numeric(max[1])]
    match.table[max[1], ] <- 0
    match.table[, max[2]] <- 0
  }

  # 05. Match all additional categories from the estimate that are
  # not in the true
  fix <- which(colnames(matches) == 0)
  new <- seq(1:NCOL(matches))
  colnames(matches)[fix] <- new[!new %in% colnames(matches)]

  # 06. Create a new column of the corrected labels and return the
  # object.
  data$estgroup <- factor(data$est, labels = colnames(matches))
  return(data)

}
