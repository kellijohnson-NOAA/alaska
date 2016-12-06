#' Reshape the results for plotting purposes
#'
#' @param data
#'
reshape_results <- function(data,
  parlabels = c("range", "rho", "sigma[e]", "sigma[epsilon]", "sigma[omega]")) {

  check <- c("replicate", "n.groups", "par", "id", "percentinc", "n.years",
    "n.groups.est")

  aa <- reshape(data, direction = "long", timevar = "par",
    varying = grep("em\\.", colnames(data)),
    drop = grep("om\\.|re\\.", colnames(data)))
  bb <- reshape(data, direction = "long", timevar = "par",
    varying = grep("re\\.", colnames(data)),
    drop = grep("om\\.|em\\.", colnames(data)))
  a <- merge(aa, bb, by = check[check %in% colnames(aa)], all = TRUE)

  aa <- reshape(data, direction = "long", timevar = "par",
      varying = grep("om\\.", colnames(data)),
      drop = grep("em\\.|re\\.", colnames(data)))
  a <- merge(a, aa, by = check[check %in% colnames(a)], all = TRUE)

  check <- length(unique(a$par)) != length(parlabels)
  if (check) stop("parlabels does not include all parameter names")
  a$par <- factor(a$par, labels = parlabels)

  a <- droplevels(a)

  return(a)
}
