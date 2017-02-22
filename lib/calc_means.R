#' Calculate the mean values for each parameter
#'
#' @param data
#' @param digits The number of significant figures you want in the
#' results for mean and MARE.
#' @param g A ggplot object with facets, of which the ranges will
#' be calculated for plotting purposes. The default is \code{NULL},
#' which allows the function to be ran with no plot.
#' @param yscalar Amount to shift the values on the y axis.
#'
calc_means <- function(data, digits = 2, ggplot = NULL,
  yscalar = 0.80) {
  MARE <- function(x) {
    median(abs(x))
  }

  a <- data[, colnames(data) %in%
    c("n_years", "par", "percentinc", "variable", "model", "om")]
  b <- data[, colnames(data) %in%
    c("n_years", "par", "percentinc", "variable", "model", "re")]
  c <- data[, colnames(data) %in%
    c("n_years", "par", "percentinc", "variable", "model", "re")]

  a <- aggregate(om ~ ., data = a, mean)
  colnames(a)[which(colnames(a) == "om")] <- "mean"

  b <- aggregate(re ~ ., data = b, MARE)
  colnames(b)[which(colnames(b) == "re")] <- "MARE"

  c <- aggregate(re ~ ., data = c, median)
  colnames(c)[which(colnames(c) == "re")] <- "MRE"

  means <- merge(a, b)
  means <- merge(means, c)

  means$MARE <- format(round(means$MARE, digits), nsmall = digits)
  means$MRE <- format(round(means$MRE, digits), nsmall = digits)

  if (!is.null(ggplot)) {
    info <- do.call(rbind, lapply(ggplot_build(ggplot)[["panel"]]$ranges,
      function(range) c(range$x.range, range$y.range)))
    info <- data.frame(info)
    info$x <- apply(info[, 1:2], 1, mean)
    info$y <- info[, 4] * yscalar
    info <- cbind(info, ggplot_build(ggplot)$panel$layout)
    means <- merge(means, info)
  }

  return(means)
}
