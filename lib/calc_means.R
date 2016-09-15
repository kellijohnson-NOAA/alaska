#' Calculate the mean values for each parameter
#'
#' @param data
#' @param digits The number of significant figures you want in the
#' results for mean and MARE.
#' @param g A ggplot object with facets, of which the ranges will
#' be calculated for plotting purposes. The default is \code{NULL},
#' which allows the function to be ran with no plot.
#'
calc_means <- function(data, digits = 2, ggplot = NULL) {
  MARE <- function(x) {
    median(abs(x))
  }
  means <- data.frame(
    "par" = unique(data$par),
    "mean" = tapply(data$om, data$par, mean),
    "MARE" = tapply(data$re, data$par, MARE))

  means$MARE <- format(round(means$MARE, digits), nsmall = digits)

  if (!is.null(ggplot)) {
    info <- do.call(rbind, lapply(ggplot_build(ggplot)[["panel"]]$ranges,
      function(range) c(range$x.range, range$y.range)))
    means$x <- apply(info[, 1:2], 1, mean)
    means$y <- info[, 4] * 0.80
  }
  return(means)
}
