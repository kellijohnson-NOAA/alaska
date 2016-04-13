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
#### Create the mesh using the projected data so all measurements are in m
###############################################################################
  ## Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  prdomain <- inla.nonconvex.hull(coordinates(data.all),
                                  -0.05,
                                  resolution = c(40, 15))

  mesh <- inla.mesh.2d(loc = coordinates(data.spp),
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
#### Set up the data
###############################################################################
  #Initial values are calculated in the model and placed in year[0]
  #year[1] is the first year of the data
  #True years are converted to 1:...
    year <- factor(data.spp$YEAR, levels = desired.years)
    levels(year) <- 1:length(desired.years)

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
    station_map <- mapval[match(station, mapval$s.unique), "map"]
  #Obtain the lat / lon (in UTM) coordinates for the used stations
    x_stations <- mesh$loc[A.locator.unique, 1]
    y_stations <- mesh$loc[A.locator.unique, 2]

###############################################################################
#### TMB
###############################################################################
  dyn.load(dynlib(file.path(dir.data, my.tmb)))
  newtonOption(smartsearch = TRUE)

    data.spatial = list(
      Y = as.vector(data.spp@data$WTCPUE),
      year = year,
      station_map = station_map,
      station_unique = station_unique,
      n_stations = length(station_unique),
      n_years = as.integer(length(desired.years) + 1), #Add a year to accommodate init pred
      G0 = spde$param.inla$M0,
      G1 = spde$param.inla$M1,
      G2 = spde$param.inla$M2
      )
    parameters.spatial = list(
      alpha = c(0.0),
      phi = 0.0,
      log_tau_E = 0.0,
      log_tau_O = 0.0,
      log_kappa = 0.0,
      log_sigma = 0.0,
      rho = 0.5,
      Epsilon_input = matrix(rnorm(spde$n.spde * (length(desired.years) + 1)),
                             nrow = spde$n.spde,
                             ncol = (length(desired.years) + 1)), #nodes x year matrix
      Omega_input = rnorm(spde$n.spde) #Omega is a vector
      )

    obj <- MakeADFun(data = data.spatial,
                     parameters = parameters.spatial,
                     random = c("Omega_input", "Epsilon_input"))

    obj$fn(obj$par)
    obj$control <- list(trace = 1, parscale = rep(1, 13), REPORT = 1,
                        reltol = 1e-12, maxit = 300)
    obj$hessian <- F
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = c(rep(-20, 2), rep(-10, 6)), #lower par bounds
                  upper = c(rep(20, 2), rep(10, 6)),   #upper par bounds
                  control = list(eval.max = 1e4, iter.max = 1e4))
    report <- try(sdreport(obj))

###############################################################################
#### Create summaries
###############################################################################

    Report_spatial <- obj$report()
    if(!("condition" %in% names(attributes(report)))) {
      opt[["summary"]] <- summary(report)
    }
    # spatial indices
    B_mean_spatial <- opt$summary[
      which(rownames(opt$summary) == "mean_abundance"), "Estimate"
      ]
    B_conf_spatial <- exp(opt$summary[
      which(rownames(opt$summary) == "log(mean_abundance)"), "Estimate"] %o%
      rep(1, 2) + opt$summary[
      which(rownames(opt$summary) == "log(mean_abundance)"), "Std. Error"] %o%
      qnorm(c(0.1, 0.9)))

  # Range of correlation (Lindgren and Rue 2013, immediately before Eq. 4)
    gmrf_range <- sqrt(8 * 1) / exp(opt$par["log_kappa"])
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_E <- 1 / sqrt(4 * pi * exp(2 * opt$par["log_tau_E"]) *
      exp(2 * opt$par["log_kappa"]))
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    marg_sd_O <- 1 / sqrt(4 * pi * exp(2 * opt$par["log_tau_O"]) *
      exp(2 * opt$par["log_kappa"]))

###############################################################################
#### rpart
###############################################################################
#Rename the stations longitude so that the output is formatted for the article
#Prune the tree according to significant differences between the relative error
#Go to http://www.statmethods.net/advstats/cart.html for more info.
longitude <- x_stations
stock <- rpart(Report_spatial$Omega ~ longitude)
stock.orig <- stock
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

###############################################################################
#### save output
###############################################################################
logfn <- create_logfile(projectName =
  paste0("alaska_", paste(tolower(substring(unlist(strsplit(desired.spp, " ")), 1, 1)),
         collapse = "")),
  projectDir = file.path(my.base, "log")
  )

  capture.output(opt, file = file.path("log", logfn), append = TRUE)
  saved <- list(
    "opt" = opt, "obj" = obj, "B_mean_spatial" = B_mean_spatial,
    "B_conf_spatial" = B_conf_spatial, "report" = report,
    "Report_spatial" = Report_spatial, "mesh" = mesh,
    "x_stations" = x_stations, "y_stations" = y_stations,
    "gmrf_range" = gmrf_range, "marg_sd_E" = marg_sd_E, "marg_sd_O" = marg_sd_O,
    "stock" = stock, "stock_all" = stock.all, "stock_orig" = stock.orig)
  save(saved, file =
    file.path("results", paste0(strsplit(logfn, ".", fixed = TRUE)[[1]][1],
              ".RData")))

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
