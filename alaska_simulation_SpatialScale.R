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
error <- 0.1
Reports <- list()
true <- c(seq(200, 2000, by = 100))
set.seed(11)
locations <- coordinates(data.all[sample(NROW(data.all), 500), ])

###############################################################################
#### Start simulation loop
###############################################################################
for (i_scale in true) {

###############################################################################
# Simulate data
###############################################################################
data_i <- Sim_Gompertz_Fn(
  n_years = diff(range(unique(data.all@data$YEAR))),
  phi = 0.0,
  SpatialScale = i_scale,
  SD_O = error,
  SD_E = error,
  SD_obs, error,
  rho = 0.5,
  logMeanDens = 1.0,
  Loc = locations,
  projection = akCRS,
  seed = 10)

meshlist <- calc_mesh(data_i[, c("Longitude", "Latitude")], NULL,
  type = "basic")

# Plot one of the true correlations from the OM
if (i_scale == true[which(true > median(true))[1]]){
  png(file.path(dir.results, paste0("SpatialScale_True_", i_scale, ".png")))
  correlations <- with(data_i[data_i$Year < 2, ],
    ncf::correlog(Longitude, Latitude, Omega,
    increment = 10, resamp = 1000, quiet = TRUE))
  plot(correlations$mean.of.class, correlations$correlation,
    ylab = "correlation", xlab = "distance (km)", las = 1)
  abline(v = i_scale, h = 0, col = "red", lty = 2)
  dev.off()
  rm(correlations)
}

###############################################################################
#### Build inputs
###############################################################################
obj <- calc_adfun(data = data_i, mesh = meshlist$mesh, tmb = my.tmb,
  variable = "zeroinflatedlnorm")

# Run optimizer
optimizer <- nlminb(obj$par, objective = obj$fn,
  gradient = obj$gr,
  lower = c(rep(-200, 5), -0.999, rep(-200, 2)),
  upper = c(rep( 200, 5), 0.999, rep(200, 2)),
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))
Reports[[length(Reports) + 1]] <- Report <- obj$report()

# Save the report
save(Report, file = file.path(dir.test,
  paste0("sim_SpatialScale_report_", i_scale, ".RData")))

rm(obj, optimizer)

} # End of i_scale loop

png(file.path(dir.results, "simulation_range.png"), units = "in",
  width = my.width[2], height = my.width[2], res = my.resolution)
temp <- data.frame("scale" = true, "range" = sapply(Reports, "[[", "Range"))
with(temp, plot(scale, range, las = 1,
  xlab = paste0("spatial scale (km) supplied to \"RMgauss(var = ", error, ")\""), ylab = ""))
temp <- lm(range ~ scale - 1, data = temp)
abline(temp)
abline(1, 1, lty = 2)
legend("topleft", legend = "1:1", lty = 2, bty = "n")
mtext(side = 2, line = 2.5, expression(range == ~ hat(sqrt(8)/kappa)))
legend("top", bty = "n",
  legend = paste("slope =", paste(round(summary(temp)$coefficients[1, 1:2], 5),
  collapse = "\nse = ")))
dev.off()

###############################################################################
#### End the file
###############################################################################

rm(data_i, error, locations, meshlist, Report, Reports, temp, true)

setwd(my.base)
# End of file
