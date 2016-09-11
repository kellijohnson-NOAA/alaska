#' Calculate input priors for the analysis
#'
#' @param input
#'
calc_priors <- function(input) {
  parameters <- list()

  parameters$alpha <- c(0.0)
  parameters$phi <- 0.0
  parameters$log_tau_E <- 1.0
  parameters$log_tau_O <- 1.0
  parameters$log_kappa <- 0.0
  parameters$rho <- 0.5
  parameters$theta_z <- c(0.0, 0.0)
  parameters$Epsilon_input <- matrix(
    rnorm(input$n_x * input$n_t), nrow = input$n_x, ncol = input$n_t)
  parameters$Omega_input <- rnorm(input$n_x)

  return(parameters)
}
