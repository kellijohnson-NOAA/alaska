#' Parse out the report into single data frame
#'
#' @param report
#'
read_report <- function(report) {

  rows <- length(report$Omega_x)
  cols <- NCOL(report$Epsilon_xt)
  lens <- length(report$t_i)

  out <- data.frame(
    # NROWS == # points in mesh
    "Site" = 1:rows,
    # scalars
    "nll" = report$jnll,
    "sd_O" = report$SigmaO,
    "sd_E" = report$SigmaE,
    "sd_obs" = report$theta_z[1],
    "theta" = report$theta_z[2],
    "phi" = report$phi,
    "alpha" = report$alpha,
    "rho" = report$rho,
    "lntau_E" = report$log_tau_E,
    "lntau_O" = report$log_tau_O,
    "SpatialScale" = report$Range,
    "lnkappa" = report$log_kappa,
    # vectors
    "eta_x" = report$eta_x,
    "Equil_x" = report$Equil_x,
    "Omega" = report$Omega_x
  )

  long <- data.frame(
    # Manipulate two matrices into long data frames
    "Year" = rep(1:cols, times = rows),
    "Site" = rep(1:rows, each = cols),
    "Epsilon" = c(report$Epsilon_xt)
  )

  some <- data.frame(
    "x_s" = report$x_s + 1,
    "t_i" = report$t_i + 1,
    "c_i" = report$c_i,
    "jnll_i" = report$jnll_i,
    "zeroinflatedlnorm" = report$log_chat_i
  )

  all <- merge(long, some, all = TRUE,
    by.x = c("Site", "Year"), by.y = c("x_s", "t_i"))
  all <- merge(all, out, all = TRUE, by = "Site")
  all <- all[order(all$Site, all$Year), ]

  return(all)
}
