###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2014-10-13
## Comments:   Polygons for the seven areas in Alaska are based on a shape file
##             provided by Annie Greig of NOAA (email: angie.greig@noaa.gov).
##             Shape file were in arcgis format.
##             The function create_areas assumes the shapefile is pocentered 
##             or projected (i.e. not in latlon format)
###############################################################################
###############################################################################

###############################################################################
#### Create the mesh
###############################################################################
  data.formesh <- eval(as.name(my.data.name))
  ## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  prdomain <- inla.nonconvex.hull(coordinates(data.formesh),
                                  -0.05,
                                  resolution = c(40, 15))

  mesh <- inla.mesh.2d(loc = coordinates(data.formesh),
                       offset = c(-0.06, -0.03),
                       cutoff = 100,
                       max.edge = c(100, 2000), 
                       boundary = prdomain)

  spde <- inla.spde2.matern(mesh, alpha = 2)

###############################################################################
#### Compile tmb file
###############################################################################
  file.remove(paste0(file.path(dir.data, my.tmb), c(".o", ".dll")))
  compile(file.path(dir.data, paste0(my.tmb, ".cpp")))

###############################################################################
#### Individual species analysis
###############################################################################

for(q in seq_along(desired.spp)){
# for(q in my.q){
  curr.spp <- desired.spp[q]
  logfile.spp <- paste(tolower(substring(unlist(strsplit(curr.spp, " ")), 
                               1, 1)), collapse = "")
  data <- eval(as.name(my.data.name))
  data <- subset(data, SID == race.num[q, 1])
  # create log file
  logfile.name <- create_logfile(projectName = paste0("alaska_", logfile.spp),
                                 projectDir = file.path(my.base, "log"))
  mesh.yes <- which(data.formesh$SID == race.num[q, 1])
  Y <- as.vector(data@data$WTCPUE)
    n_years <- length(desired.years)
    year_0 <- year <- factor(data$YEAR, levels = desired.years)
    levels(year) <- 1:n_years
    levels(year_0) <- 1:n_years - 1

  A.locator <- mesh$idx$loc[mesh.yes]
    A.locator.unique <- unique(A.locator)[order(unique(A.locator))] 
    station <- A.locator - 1
    station_unique  <- A.locator.unique - 1
    n_stations <- length(station_unique)
      
      newlevels <- 1:spde$n.spde
      station_map <- factor(A.locator, levels = newlevels)
      names(newlevels) <- NA
      names(newlevels)[A.locator.unique] <- 1:(n_stations)
      station_map <- as.numeric(names(newlevels)[station_map]) - 1
    x_stations <- mesh$loc[A.locator.unique, 1]
    y_stations <- mesh$loc[A.locator.unique, 2]


###############################################################################
#### TMB
###############################################################################
  dyn.load(dynlib(file.path(dir.data, my.tmb)))
  newtonOption(smartsearch = TRUE)

if(my.tmb == "gompertz_kfj"){
    data.spatial = list(Y = Y, 
                        year = year,
                        station_map = station_map,
                        station_unique = station_unique,
                        n_stations = n_stations,
                        n_years = as.integer(n_years + 1),
                        G0 = spde$param.inla$M0, 
                        G1 = spde$param.inla$M1, 
                        G2 = spde$param.inla$M2)
    parameters.spatial = list(alpha = c(0.0),
                              phi = 0.0,
                              log_tau_E = 0.0,
                              log_tau_O = 0.0,
                              log_kappa = 0.0,
                              log_sigma = 0.0,
                              rho = 0.5, 
                              Epsilon_input = matrix(rnorm(spde$n.spde * (n_years + 1)),
                                                     nrow = spde$n.spde,
                                                     ncol = (n_years + 1)), 
                              Omega_input = rnorm(spde$n.spde))

    obj <- MakeADFun(data = data.spatial, 
                     parameters = parameters.spatial,
                     random = c("Omega_input", "Epsilon_input"))
}

    obj$fn(obj$par)
    obj$control <- list(trace = 1, parscale = rep(1, 13), REPORT = 1,
                                reltol = 1e-12, maxit = 300)
    obj$hessian <- F
    opt <- nlminb(obj$par, obj$fn, obj$gr, 
                  lower = c(rep(-20, 2), rep(-10, 6)), 
                  upper = c(rep(20, 2), rep(10, 6)),
                  control = list(eval.max = 1e4, iter.max = 1e4))
    report <- try(sdreport(obj))
    Report_spatial <- obj$report()
    if(!("condition" %in% names(attributes(report)))) {
      opt[["summary"]] <- summary(report)
    } 
    # spatial indices
    B_mean_spatial <- opt$summary[which(rownames(opt$summary) == "mean_abundance"),
                                  "Estimate"]
    B_conf_spatial <- exp(opt$summary[which(rownames(opt$summary) == "log(mean_abundance)"),
                                      "Estimate"] %o% rep(1, 2) + 
                      opt$summary[which(rownames(opt$summary) == "log(mean_abundance)"),
                                  "Std. Error"] %o% qnorm(c(0.1, 0.9)))

  # Range of correlation (Lindgren and Rue 2013, immediately before Eq. 4)
    gmrf_range <- sqrt(8*1)/exp(opt$par["log_kappa"])
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_E <- 1 / sqrt(4*pi*exp(2*opt$par["log_tau_E"])*exp(2*opt$par["log_kappa"]))
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_O <- 1 / sqrt(4*pi*exp(2*opt$par["log_tau_O"])*exp(2*opt$par["log_kappa"]))

###############################################################################
#### rpart
###############################################################################
longitude <- x_stations
stock <- rpart(Report_spatial$Omega ~ longitude)
stock.original <- stock
splits <- rep(NA, (dim(stock$cptable)[1] - 1))

if(length(splits) > 1){
  for(i in 1:(dim(stock$cptable)[1] - 1)){
    x1 <- stock$cptable[i, "xerror"]
    x2 <- stock$cptable[i + 1, "xerror"]
    meandiff <- x1 - x2
    denom <- sqrt(stock$cptable[i, "xstd"]^2 + 
                  stock$cptable[i + 1, "xstd"]^2)
    splits[i] <- ifelse(0 < meandiff + 1.96 * denom & 
                        0 > meandiff - 1.96 * denom, 
                        "same", "different")
}
  cp.choice <- stock$cptable[max(which(splits == "different")) + 1, "CP"]
  temp.coords <- data.frame(x = stock$splits[, "index"], 
                            y = rep(mean(y_stations), dim(stock$splits)[1]))
  coordinates(temp.coords) <- ~ x + y
  proj4string(temp.coords) <- akCRS
  temp.coords <- spTransform(temp.coords, llCRS)
  stock$splits[, "index"] <- temp.coords@coords[, "x"]
  rownames(stock$splits) <- rep("longitude", dim(stock$splits)[1])
  stock.all <- stock
  stock <- prune(stock, cp = cp.choice)

}


## Go to http://www.statmethods.net/advstats/cart.html for more info.
###############################################################################
#### save output
###############################################################################
  capture.output(opt, file = file.path("log", logfile.name), append = TRUE)
  saved <- list("opt" = opt, "obj" = obj,
                "B_mean_spatial" = B_mean_spatial, "B_conf_spatial" = B_conf_spatial, 
                "report" = report, "Report_spatial" = Report_spatial, 
                "mesh" = mesh, "x_stations" = x_stations, "y_stations" = y_stations,
                "gmrf_range" = gmrf_range, "marg_sd_E" = marg_sd_E,
                "marg_sd_O" = marg_sd_O, 
                "stock" = stock, "stock_all" = stock.all, "stock_orig" = stock.orig)
  save(saved, file = file.path("results", paste0(strsplit(logfile.name, ".", 
                                                          fixed = TRUE)[[1]][1],
                                                 ".RData")))


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
}

#EndOfFile