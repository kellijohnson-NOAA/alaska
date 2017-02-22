#' Calculate the relative error of a parameter
#'
#' @details The function relies on the column names to be matched
#' for the operating and estimation methods,
#' i.e., \code{"Omega_om"} and \code{"Omega_em"}.
#'
#' @param data
#' @param return A character value specifying what to return.
#' If \code{return = "all"} the orginal data frame is bound to the
#' relative error calculations, and if \code{return = "re"} then only
#' the newly determined measurements of relative error are returned.
#'
#' @author Kelli Faye Johnson
#'
#' @examples
#' temp <- data.frame("Omega_om" = 2:3, "Omega_em" = c(3,3))
#' calc_RE(temp)
calc_RE <- function(data, return = c("all", "re")) {
  # Argument checking
  return <- match.arg(return, c("all", "re"), several.ok = FALSE)
  if (is.vector(data)) {
    stop("calc_RE only works on data frames.")
  }

  # Get matching om and em names
  names <- colnames(data)
  # Change _ to __ for om and em
  names <- gsub("_(?=[[:alpha:]]m)", "__", names, perl = TRUE)

  # Find the matching column
  match <- data.frame(do.call("rbind", strsplit(names, "__")))
  match$match <- row.names(match)
  for (ii in 1:NROW(match)) {
    if (match[ii, 2] %in% c("om", "em")) {
      temp <- which(match[, 1] == match[ii, 1])
      match$match[ii] <- temp[!temp %in% ii]
    } else { next }
  }
  match <- droplevels(match[match[, 2] %in% c("om", "em"), ])

  # Determine if an om par has a value of 0
  oms <- data[, grep("[[:alpha:]]_om", colnames(data))]
  oms <- colnames(oms)[which(apply(oms, 2, mean, na.rm = TRUE) == 0)]
  oms <- grep(paste(gsub("_om", "", oms), collapse = "|"),
    colnames(data))
  for (ii in oms) {
    data[, ii] <- exp(data[, ii])
  }

  # Calculate the relative error (RE)
  re <- tapply(match$match, match[, 1], function(x) {
    x <- as.numeric(x)
    (data[, x[2]] - data[, x[1]]) / data[, x[2]]
   })
  if (!is.list(re)) {
    re <- data.frame(t(re))
  } else {
    re <- t(do.call("rbind", re))
  }
  colnames(re) <- paste0("re_", colnames(re))

  if (return == "re") return(re)
  if (return == "all") return(cbind(data, re))

}
