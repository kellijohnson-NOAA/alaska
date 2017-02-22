#' Reshape the results for plotting purposes
#'
#' @param data
#'
reshape_results <- function(data,
  parlabels = c("alpha", "epsilon", "ln(size)",
    "omega", "phi", "rho", "sigma[epsilon]",
    "sigma[omega]", "sigma[obs]", "scale", "zeroinflatedlnorm")) {

  check <- c("replicate", "subpopulations", "percentinc", "n_years",
    "par", "id", "model", "variable")
  if (!all(check[1:4] %in% colnames(data))) {
    stop("not all of the fixed columns are present in data")
  }
  keep <- c(which(colnames(data) %in% check),
    grep("_om|_em|re_", colnames(data)))
  data <- data[, keep]

  aa <- reshape(data, direction = "long", timevar = "par",
    v.names = "em",
    varying = list(grep("_em", colnames(data))),
    times = gsub("_em", "", colnames(data)[grep("_em", colnames(data))]),
    drop = grep("_om|re_", colnames(data)))
  bb <- reshape(data, direction = "long", timevar = "par",
    v.names = "om",
    varying = list(grep("_om", colnames(data))),
    times = gsub("_om", "", colnames(data)[grep("_om", colnames(data))]),
    drop = grep("_em|re_", colnames(data)))
  a <- merge(aa, bb, by = check[check %in% colnames(aa)], all = TRUE)

  if (!any(grepl("re_", colnames(data)))) {
    stop("Must calculate relative error before reshaping the data.")
  }
  aa <- reshape(data, direction = "long", timevar = "par",
    v.names = "re",
      varying = list(grep("re_", colnames(data))),
      times = gsub("re_", "", colnames(data)[grep("re_", colnames(data))]),
      drop = grep("_om|_em", colnames(data)))
  a <- merge(a, aa, by = check[check %in% colnames(a)], all = TRUE)

  check <- length(unique(a$par)) != length(parlabels)
  if (check) stop("parlabels does not include all parameter names")
  a$par <- factor(a$par, labels = parlabels)

  a <- droplevels(a)

  return(a)
}
