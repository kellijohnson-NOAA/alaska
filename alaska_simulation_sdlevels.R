###############################################################################
###############################################################################
## Purpose:    Code for Appendix figure with different levels of process and
##             omega standard deviations for spatial variation.
##
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-10-31
## Comments:
###############################################################################
###############################################################################

###############################################################################
# Plot what changing the levels of process and observation error does to the
# dynamics of the population
###############################################################################
errors <- c("level" = c(0.01, 0.50))
pdata <- NULL
for (omega in seq_along(errors)) {
for (process in seq_along(errors)) {
  pdata <- rbind(pdata, data.frame(Sim_Gompertz_Fn(
    n_years = diff(range(unique(data.all@data$YEAR))),
    SpatialScale = 0.25,
    SD_O = errors[omega],
    SD_E = errors[process],
    SD_obs = 0.05,
    logMeanDens = 0.2,
    rho = 0.5,
    phi = 0.0,
    n_stations = 200,
    Loc = NULL,
    projection = akCRS,
    seed = 10),
    "process" = paste("sigma[epsilon] ==", errors[process]),
    "omega" = paste("sigma[omega] ==", errors[omega])))
  if (is.null(pdata)) pdata <- pdata$DF
  else pdata <- rbind(pdata, pdata$DF)
} # End cols
} # End rows

png(file.path(dir.results, "simulation_sdlevels.png"),
  width = my.width[2], height = my.width[2],
  units = "in", res = my.resolution,
  )

pdata <- pdata[pdata$Year == max(pdata$Year), ]
colnames(pdata)[which(colnames(pdata) == "Simulated_example")] <-
 "count"
ggplot(data = pdata, aes(y = Latitude, x = Longitude)) +
facet_grid(process ~ omega, labeller = label_parsed) +
geom_point(aes(size = count), fill = NA, shape = 21) +
scale_size(limits = c(1, max(pdata$count))) +
xlab("easting") + ylab("northing") +
coord_equal() + my.theme +
theme(
  legend.key = element_rect(fill = NA, colour = NA, size = 0.25),
  panel.border = element_rect(colour = "black", fill = NA, size = 0.1),
  legend.title = element_text(size = 15),
  strip.text.y = element_text(size = 15),
  strip.text.x = element_text(size = 15))

dev.off()
rm(errors, pdata, omega, process)
