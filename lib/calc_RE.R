#' Calculate the relative error of a parameter
#'
#' @param data
#'
calc_RE <- function(data) {
  # Get matching om and em names
  if (is.vector(data)) {
    data <- data.frame(t(data))
  }
  names <- colnames(data)

  # Find the matching row for each OM
  match <- data.frame(do.call("rbind", strsplit(names, "_")))
  match$match <- row.names(match)
  for (ii in 1:NROW(match)) {
    if (match[ii, 1] == "om") {
      temp <- which(match[, 2] == match[ii, 2])
      match$match[ii] <- ifelse(length(temp) > 1, tail(temp, 1), NA)
    } else { if (match[ii, 1] == "em") {
      match$match[ii] <- which(match[, 2] == match[ii, 2])[1]
      } else { next }}
  }
  match <- droplevels(match[match[, 1] %in% c("om", "em"), ])

  # Calculate the relative error (RE)
  re <- tapply(match$match, match[, 2], function(x) {
    x <- as.numeric(x)
    (data[, x[2]] - data[, x[1]]) / data[, x[2]]
   })
  if (!is.list(re)) {
    re <- data.frame(t(re))
  } else {
    re <- t(do.call("rbind", re))
  }
  colnames(re) <- paste0("re_", colnames(re))

  return(re)
}
