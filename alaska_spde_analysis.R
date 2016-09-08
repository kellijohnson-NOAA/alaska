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
#### Compile tmb file
###############################################################################
file.remove(paste0(file.path(dir.data, my.tmb), c(".o", ".dll")))
compile(file.path(dir.data, paste0(my.tmb, ".cpp")))
dyn.load(dynlib(file.path(dir.data, my.tmb)))

###############################################################################
#### Set up the data
###############################################################################
  #Initial values are calculated in the model and placed in t_i[0]
  #t_i[1] is the first year of the data

  #A.locator is the location of each data point on the mesh
  #Only some of the mesh nodes are actually filled
  #Filled nodes == A.locator.unique
  #station is used to map 1:numberofnodes to 0:numberofnodes-1
  #station_map creates a vector to map the 0:numberofusednodes-1
  #to 0:numberofnodes-1 because the cpp code only predicts for used nodes
  A.locator <- mesh$idx$loc
    A.locator.unique <- unique(A.locator)[order(unique(A.locator))]
    station <- A.locator - 1
    station_unique  <- A.locator.unique - 1
    mapval <- data.frame(s.unique = station_unique,
                         map = seq(0, (length(A.locator.unique) - 1)))
  #Obtain the lat / lon (in UTM) coordinates for the used stations
    x_stations <- mesh$loc[A.locator.unique, 1]
    y_stations <- mesh$loc[A.locator.unique, 2]

###############################################################################
#### TMB
###############################################################################
# Build inputs
  X_xp <- matrix(1, ncol = 1, nrow = spde$n.spde)

  data <- list(
    Options_vec = 1,
    n_i = NROW(data.all@data),
    n_x = spde$n.spde,
    n_t = length(desired.years),
    n_p = NCOL(X_xp),

    x_s = mesh$idx$loc - 1,

    c_i = as.vector(data.all@data$WTCPUE),
    s_i = mapval[match(station, mapval$s.unique), "map"],
    t_i = factor(data.all@data$YEAR, levels = desired.years,
      labels = 1:length(desired.years) - 1),

    X_xp = X_xp,

    G0 = spde$param.inla$M0,
    G1 = spde$param.inla$M1,
    G2 = spde$param.inla$M2
    )

  parameters = list(
    alpha = c(0.0),
    phi = 0.0,
    log_tau_E = 1.0,
    log_tau_O = 1.0,
    log_kappa = 0.0,
    rho = 0.5,
    theta_z = c(0, 0),
    Epsilon_input = matrix(rnorm(spde$n.spde * data$n_t),
      nrow = spde$n.spde, ncol = data$n_t),
    Omega_input = rnorm(spde$n.spde)
    )

  obj <- MakeADFun(data = data,
    parameters = parameters,
    random = c("Epsilon_input", "Omega_input"),
    map = NULL,
    hessian = FALSE,
    DLL = my.tmb)

newtonOption(obj, smartsearch = FALSE)

  opt <- nlminb(obj$par, obj$fn, gradient = obj$gr,
    lower = c(rep(-200, 5), -0.999, rep(-200, 2)), #lower par bounds
    upper = c(rep(200, 5), 0.999, rep(200, 2)),    #upper par bounds
    control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))
  opt[["final_gradient"]] <- obj$gr(opt$par)

  # Obtain standard errors
  Report <- obj$report()
  rm(report)
  report <- try(sdreport(obj))
  unlist(Report[c('Range', 'SigmaO', 'SigmaE', 'rho', 'theta_z')])

###############################################################################
#### Clustering
###############################################################################
# Get the points of the nodes from the mesh
plotgroups <- data.frame(
  "x" = mesh$loc[, 1], "y" = mesh$loc[, 2],
  "latitude" = round(mesh$loc[, 1], 1),
  "omega" = Report[["Omega_x"]],
  "clustering" = NA)
coordinates(plotgroups) <- ~ x + y
proj4string(plotgroups) <- akCRS

# Find points in the main frame of the mesh
localboundaries <- findlocal(mesh)
plotall <- plotgroups
plotgroups <- plotgroups[localboundaries, ]

# Calculate the clusters
est <- list()
est$kdist <- dist(Report[["Omega_x"]][localboundaries],
  method = "euclidean")

# Find the optimum number of groups
est$cluster.ksearch <- cluster_validity(dist.file = est$kdist,
  tot_k = ceiling(sqrt(attr(est$kdist, "Size") - 1)),
  plot_sils = TRUE, file = file.path(dir.results, "est_ksearch.png"))
est$clustermax <- sapply(est$cluster.ksearch[-1], which.max)
est$clustermin <- sapply(est$cluster.ksearch[-1], which.min)
est$cluster.use <- cluster::pam(x = est$kdist,
  k = which.max(est$cluster.ksearch$Hubert.gamma) + 1)
# Assign the clusters to the spatial points data frame
plotgroups$clustering <- est$cluster.use$clustering

# Estimate the 2D groupings from the clustering
est$spatial <- SPODT::spodt(omega ~ 1,
  data = plotgroups, level.max = 3)
est$lines <- SPODT::spodtSpatialLines(est$spatial,
  data = plotgroups)

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
