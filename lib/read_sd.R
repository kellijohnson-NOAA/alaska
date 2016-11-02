#' Create a single data frame of uncertainty measures
#'
#' @param sd
#'
read_sd <- function(sd) {

  out <- data.frame(
    "par" = c(names(sd$value), names(sd$par.fixed)),
    "est" = c(sd$value, sd$par.fixed),
    "var" = c(diag(sd$cov), diag(sd$cov.fixed)))
  return(out)
}
