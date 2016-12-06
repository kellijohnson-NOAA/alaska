#' Read the parameter values from the results file
#'
#' @param results
#'
read_pars <- function(results) {
  pars <- c(
    "replicate" = results$DF$replicate[1],
    "subpopulations" = results$DF$subpopulations[1],
    "percentinc" = results$DF$percentinc[1],
    "n_years" = length(unique(results$DF$Year)),
    "om_Range" = results$DF$SpatialScale[1],
    "om_rho" = 0.5,
    "om_SigmaO" = results$DF$sd_O[1],
    "om_SigmaE" = results$DF$sd_E[1],
    "om_Sigmae" = results$DF$sd_obs[1],
    "em_Range" = results$Range,
    "em_rho" = results$rho,
    "em_SigmaO" = results$SigmaO,
    "em_SigmaE" = results$SigmaE,
    "em_Sigmae" = exp(results$theta_z[1])
  )

  return(pars)
}
