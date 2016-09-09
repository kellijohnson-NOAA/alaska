#' Calculate the relative error of a parameter
#'
#' @param data
#'
calc_RE <- function(data) {
  # Get matching om and em names
  if (is.vector(data)) {
    data <- data.frame(t(pars))
  }
  names <- colnames(data)

  # Find the matching row for each OM
  match <- data.frame(do.call("rbind", strsplit(names, "_")))
  match$match <- rep(NA, NROW(match))
  for (ii in 1:NROW(match)) {
    if (match[ii, 1] == "om") {
      temp <- which(match[, 2] == match[ii, 2])
      match$match[ii] <- ifelse(length(temp) > 1, tail(temp, 1), NA)
    } else { next }
  }

  # Calculate the relative error (RE)
  re <- tapply(match[,3], match[,2], function(x) {
    (pars[, x[1]] - pars[, x[2]]) / pars[, x[1]]
  })
  if (!is.list(re) == 1) {
    re <- data.frame(t(re))
  } else {
    re <- do.call("rbind", re)
  }
  colnames(re) <- paste0("re_", colnames(re))

  return(re)
}
