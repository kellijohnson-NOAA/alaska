###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2015-01-05
## Comments:
###############################################################################
###############################################################################

###############################################################################
#### Set up the data
###############################################################################
  #Initial values are calculated in the model and placed in t_i[0]
  #t_i[1] is the first year of the data


  #Obtain the lat / lon (in UTM) coordinates for the used stations
    x_stations <- mesh$loc[unique(mesh$idx$loc)[order(unique(mesh$idx$loc))], 1]
    y_stations <- mesh$loc[unique(mesh$idx$loc)[order(unique(mesh$idx$loc))], 2]

###############################################################################
#### TMB
###############################################################################
# Build inputs
  data <- list(
    Options_vec = 1,
    n_i = NROW(data.all@data),
    n_x = spde$n.spde,
    n_t = length(desired.years),
    n_p = 1,

    x_s = mesh$idx$loc - 1,

    c_i = as.vector(data.all@data$WTCPUE),
    s_i = calc_meshmap(mesh),
    t_i = factor(data.all@data$YEAR, levels = desired.years,
      labels = 1:length(desired.years) - 1),

    X_xp = matrix(1, ncol = 1, nrow = spde$n.spde),

    G0 = spde$param.inla$M0,
    G1 = spde$param.inla$M1,
    G2 = spde$param.inla$M2
    )

  parameters <- calc_priors(data)

  obj <- MakeADFun(data = data,
    parameters = parameters,
    random = c("Epsilon_input", "Omega_input"),
    map = NULL, # map is always NULL b/c using rates instead of count
    hessian = FALSE,
    DLL = my.tmb,
    silent = TRUE)

newtonOption(obj, smartsearch = FALSE)

  opt <- nlminb(obj$par, obj$fn, gradient = obj$gr,
    lower = c(rep(-200, 5), -0.999, rep(-200, 2)), #lower par bounds
    upper = c(rep(200, 5), 0.999, rep(200, 2)),    #upper par bounds
    control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))
  opt[["final_gradient"]] <- obj$gr(opt$par)

  # Obtain standard errors
  Report <- obj$report()
  report <- try(sdreport(obj))
  unlist(Report[c('Range', 'SigmaO', 'SigmaE', 'rho', 'theta_z')])
  save(Report, report, file = file.path(dir.results, "analysis.RData"))

###############################################################################
#### Clustering
###############################################################################
# Get the points of the nodes from the mesh
plotgroups <- data.frame(
  "x" = mesh$loc[, 1], "y" = mesh$loc[, 2],
  "omega" = Report[["Omega_x"]],
  "est" = NA)
coordinates(plotgroups) <- ~ x + y
proj4string(plotgroups) <- akCRS
plotall <- plotgroups

# Find points in the main frame of the mesh
plotgroups <- plotgroups[findlocal(mesh), ]

# Calculate the clusters
    export_prm(file_in = my.prm, dir = getwd(), file = "analysis",
      shape = "circle", size = 100)
    export_geo(points = spTransform(plotgroups, llCRS),
      dir = getwd(), file = "analysis", projection = llCRS)
    # Double check the direction of slashes and call system
    system(paste(
      gsub("/", "\\\\", my.SaTScan),
      file.path(gsub("/", "\\\\", getwd()), paste0("analysis", ".prm"), fsep = "\\")
    ), show.output.on.console = FALSE)
    # Read in the results
    textfile <- read_txt(dir = getwd(), file = "analysis")
    shape <- readShapePoly(file.path(getwd(), "analysis.col"))
    projection(shape) <- llCRS

# Get p values
if (is.data.frame(textfile$coordinates) | is.matrix(textfile$coordinates)) {
  pvalues <- as.numeric(sapply(strsplit(
    textfile$coordinates[grep("P-", textfile$coordinates[, 1]), ],
    ":"), "[[", 2))
}

load(file.path(dir.data, "maps.eez.RData"))

png(file.path(dir.results, "analysis_clusters.png"),
  width = my.width[2], height = my.height.map, res = my.resolution, units = "in")
plot(maps.eez)
r4kfj::llgridlines(maps.eez, recenter = TRUE, lty = 1, col = col.gridlines)
plot(spTransform(plotgroups, efhCRS), cex = 5 * (abs(plotgroups$omega)),
  pch = ifelse(plotgroups$omega >= 0, 15, 0), add = TRUE,
  col = rgb(0, 0, 0, 0.4))
for (it_ in seq_along(shape)) {
  if (pvalues[it_] > 0.05) next
  plot(spTransform(subset(shape, CLUSTER == it_), efhCRS), add = TRUE)
}
dev.off()



###############################################################################
#### Create summaries
###############################################################################

    if(!("condition" %in% names(attributes(report)))) {
      opt[["summary"]] <- summary(report)
    }

  # Range of correlation (Lindgren and Rue 2013, immediately before Eq. 4)
    gmrf_range <- sqrt(8 * 1) / exp(opt$par["log_kappa"])
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_E <- 1 / sqrt(4 * pi * exp(2 * opt$par["log_tau_E"]) *
      exp(2 * opt$par["log_kappa"]))
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_O <- 1 / sqrt(4 * pi * exp(2 * opt$par["log_tau_O"]) *
      exp(2 * opt$par["log_kappa"]))

###############################################################################
#### unload the cpp file
###############################################################################
  # First attempt
  firsttry <- try(dyn.unload(dynlib(file.path(dir.data, my.tmb))),
                  silent = TRUE)
  # Garbage Collection
  if(is(firsttry, "try-error")) gc()
  # Second attempt
  secondtry <- try(dyn.unload(dynlib(file.path(dir.data, my.tmb))),
                   silent = TRUE)
  # Verify that second attempt works
  #getLoadedDLLs()

#EndOfFile
