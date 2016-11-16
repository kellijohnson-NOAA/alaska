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
    "sd_O_est" = report$SigmaO,
    "sd_E_est" = report$SigmaE,
    "sd_obs_est" = report$theta_z[1],
    "theta_est" = report$theta_z[2],
    "phi_est" = report$phi,
    "alpha_est" = report$alpha,
    "rho_est" = report$rho,
    "lntau_E_est" = report$log_tau_E,
    "lntau_O_est" = report$log_tau_O,
    "SpatialScale_est" = report$Range,
    "lnkappa_est" = report$log_kappa,
    # vectors
    "eta_x" = report$eta_x,
    "Equil_x" = report$Equil_x,
    "Omega_est" = report$Omega_x
  )

  long <- data.frame(
    # Manipulate two matrices into long data frames
    "Year" = rep(1:cols, times = rows),
    "Site" = rep(1:rows, each = cols),
    "Epsilon_est" = c(report$Epsilon_xt)
  )

  some <- data.frame(
    "x_s" = rep(report$x_s + 1, each = cols),
    "t_i" = report$t_i + 1,
    "c_i" = report$c_i,
    "jnll_i" = report$jnll_i,
    "log_chat_i" = report$log_chat_i
  )

  all <- merge(long, some, all = TRUE,
    by.x = c("Site", "Year"), by.y = c("x_s", "t_i"))
  all <- merge(all, out, all = TRUE, by = "Site")
  all <- all[order(all$Site, all$Year), ]

  return(all)
}
