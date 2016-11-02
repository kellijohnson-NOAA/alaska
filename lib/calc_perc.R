#' Calculate the percent increase in alpha across a scenario.
#'
#' @param scenario A data frame with a scenario.
#'
calc_perc <- function(scenario) {
  if(length(scenario) == 1) return(0)
  if(NCOL(scenario) > 2) stop("Need to fix calc_perc",
    "to accommodate more than two groups.")
  else {return(
    apply(scenario, 1, diff) /
    apply(scenario, 1, min)[it_alpha]
    )
  }
}
