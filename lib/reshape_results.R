#' Reshape the results for plotting purposes
#'
#' @param data
#'
reshape_results <- function(data,
  parlabels = c("range", "rho", "sigma[Epsilon]", "sigma[Omega]")) {

  data <- merge(
    reshape(data, direction = "long", timevar = "par",
      varying = grep("om\\.", colnames(data)),
      drop = grep("em\\.|re\\.", colnames(data))),
    merge(
    reshape(data, direction = "long", timevar = "par",
      varying = grep("em\\.", colnames(data)),
      drop = grep("om\\.|re\\.", colnames(data))),
    reshape(data, direction = "long", timevar = "par",
      varying = grep("re\\.", colnames(data)),
      drop = grep("om\\.|em\\.", colnames(data))),
    by = c("replicate", "n.groups", "par", "id", "percentinc")
    ),
    by = c("replicate", "n.groups", "par", "id", "percentinc"))
  data$par <- factor(data$par, labels = parlabels)

  return(data)
}
