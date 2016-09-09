#' Read the parameter values from the results file
#'
#' @param results
#'
read_pars <- function(results) {
  pars <- c(
    "om_Range" = results$input$SpatialScale,
    "om_rho" = results$input$rho,
    "om_SigmaO" = results$input$SD_O,
    "om_SigmaE" = results$input$SD_O,
    "em_Range" = results$Range,
    "em_rho" = results$rho,
    "em_SigmaO" = results$SigmaO,
    "em_SigmaE" = results$SigmaE
  )

  return(pars)
}
