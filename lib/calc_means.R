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

  if ("n.years" %in% colnames(data)) {
    a <- aggregate(om ~ par + percentinc + n.years, data = data, mean)
  } else a <- aggregate(om ~ par + percentinc, data = data, mean)
  colnames(a)[which(colnames(a) == "om")] <- "mean"

  if ("n.years" %in% colnames(data)) {
    b <- aggregate(re ~ par + percentinc + n.years, data = data, MARE)
  } else b <- aggregate(re ~ par + percentinc, data = data, MARE)
  colnames(b)[which(colnames(b) == "re")] <- "MARE"

  means <- merge(a, b)

  means$MARE <- format(round(means$MARE, digits), nsmall = digits)

  if (!is.null(ggplot)) {
    info <- do.call(rbind, lapply(ggplot_build(ggplot)[["panel"]]$ranges,
      function(range) c(range$x.range, range$y.range)))
    info <- data.frame(info)
    info$par <- rep(levels(means$par),
      ifelse(length(unique(means$percentinc)) > 1,
        NROW(info) / length(levels(means$par)), 1))
    info$percentinc <- rep(unique(means$percentinc),
      each = NROW(info) / length(unique(means$percentinc)))
    info$x <- apply(info[, 1:2], 1, mean)
    info$y <- info[, 4] * yscalar
    means <- merge(means, info[, c("par", "percentinc", "x", "y")],
      by = c("par", "percentinc"))
  }

  return(means)
}
