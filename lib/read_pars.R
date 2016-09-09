#' Read the parameter values from the results file
#'
#' @param results
#'
read_pars <- function(results) {
  pars <- c(
    "om_Range" = results$SpatialScale,
    "om_rho" = results$om_rho,
    "om_SigmaO" = results$SD_O,
    "om_SigmaE" = results$SD_O,
    "em_Range" = results$Range,
    "em_rho" = results$rho,
    "em_SigmaO" = results$SigmaO,
    "em_SigmaE" = results$SigmaE,
  )

  return(pars)
}
