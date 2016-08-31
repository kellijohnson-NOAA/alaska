###############################################################################
###############################################################################
## Purpose:    Simulate data for a spatial gompertz model
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-08-30
## Comments:
###############################################################################
###############################################################################
dir.test <- file.path(my.base, "simulation_SpatialScale")
dir.create(dir.test, showWarnings = FALSE)
setwd(dir.test)

###############################################################################
# Provide simulation inputs
###############################################################################
error <- 0.01
Reports <- list()
true <- c(seq(200, 2000, by = 50))

###############################################################################
#### Unload and load cpp
###############################################################################
calc_cpp(cpp = my.tmb, loc = dir.data)

for (i_scale in true) {
  set.seed(10)

###############################################################################
# Simulate data
###############################################################################
sim_data <- Sim_Gompertz_Fn(
  n_years = length(unique(data.all@data$YEAR)),
  SpatialScale = i_scale,
  SD_O = error,
  SD_E = error,
  SD_extra = 0.0,
  rho = 0.5,
  logMeanDens = c(2.0),
  phi = 0.0,
  Loc = coordinates(data.all[sample(NROW(data.all), 1000), ]),
  projection = akCRS,
  weightvals = c(
    mean(data.all$WTCPUE / data.all$NUMCPUE, na.rm = TRUE),
    sd(data.all$WTCPUE / data.all$NUMCPUE, na.rm = TRUE)))

meshlist <- calc_mesh(sim_data$DF[, c("Longitude", "Latitude")])

# ggplot(data = cbind(sim_data$DF, "group" = factor(rep(sim_data$group, each = 22))),
#   aes(x = Year, y = Simulated_example)) +
#   geom_point(aes(col = group)) +
#   ylab("Count")

if (i_scale == 1000){
  png(file.path(dir.results, paste0("SpatialScale_True_", i_scale, ".png")))
  correlations <- ncf::correlog(sim_data[["Loc"]][,1],
    sim_data[["Loc"]][,2], sim_data[["Omega"]],
    increment = 10, resamp = 1000, quiet = TRUE)
  par(las = 1)
  plot(correlations$mean.of.class, correlations$correlation,
    ylab = "correlation", xlab = "distance (km)")
  abline(v = i_scale, h = 0, col = "red", lty = 2)
  dev.off()
}

###############################################################################
#### Build inputs
###############################################################################
Data <- list(
  Options_vec = ifelse(all(round(
    sim_data$DF[, "Simulated_example"]) == sim_data$DF[, "Simulated_example"]),
    0, 1),
  n_i = NROW(sim_data$DF),
  n_x = meshlist$mesh$n,
  n_t = max(sim_data$DF$Year),
  n_p = 1,
  x_s = meshlist$mesh$idx$loc - 1,
  c_i = sim_data$DF[, "Simulated_example"],
  s_i = sim_data$DF[, "Site"] - 1,
  t_i = sim_data$DF[, "Year"] - 1,
  X_xp = matrix(1, ncol = 1, nrow = meshlist$mesh$n),
  G0 = meshlist$spde$param.inla$M0,
  G1 = meshlist$spde$param.inla$M1,
  G2 = meshlist$spde$param.inla$M2)
Parameters <- list(
  alpha = c(0.0),
  phi = 0.0,
  log_tau_E = 1.0,
  log_tau_O = 1.0,
  log_kappa = 0.0,
  rho = 0.5,
  theta_z = c(0.0, 0.0),
  Epsilon_input = matrix(rnorm(meshlist$mesh$n * Data$n_t),
    nrow = meshlist$mesh$n, ncol = Data$n_t),
  Omega_input = rnorm(meshlist$mesh$n))
Random <- c("Epsilon_input", "Omega_input")
Map <- NULL
if (Data$Options_vec == 0) {
  Map[["theta_z"]] <- factor(c(NA, NA))
}

# Make object
obj <- MakeADFun(
  data = Data,
  parameters = Parameters,
  random = Random,
  map = Map,
  hessian = FALSE,
  DLL = my.tmb)

# Run optimizer
optimizer <- nlminb(obj$par, objective = obj$fn,
  gradient = obj$gr,
  lower = c(rep(-200, 5), -0.999, rep(-200, 2)),
  upper = c(rep( 200, 5), 0.999, rep(200, 2)),
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))
Reports[[length(Reports) + 1]] <- Report <- obj$report()

# Save the report
save(Report, sim_data, Data, Parameters, file = file.path(dir.test,
  paste0("sim_SpatialScale_report_", i_scale, ".RData")))

rm(obj)

} # End of i_scale loop

png(file.path(dir.results, "SpatialScale_Range.png"), units = "in",
  width = my.width[2], height = my.width[2], res = my.resolution)
temp <- data.frame("scale" = true, "range" = sapply(Reports, "[[", "Range"))
# Select only those models with a rho having < 40% RE
temp <- temp[which(sapply(Reports, "[[", "rho") < (0.3 * 0.5 + 0.5)), ]
with(temp, plot(scale, range, las = 1,
  xlab = "spatial scale supplied to \"RMgauss()\" (km)", ylab = ""))
temp <- lm(range ~ scale, data = temp)
abline(temp)
mtext(side = 2, line = 2.2, expression(range == ~ hat(sqrt(8)/kappa)))
legend("bottomright", bty = "n",
  legend = paste("slope =", paste(round(summary(temp)$coefficients[2, 1:2], 5),
  collapse = ", se = ")))
rm(temp)
dev.off()

setwd(my.base)
# End of file
