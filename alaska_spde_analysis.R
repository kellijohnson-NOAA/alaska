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
# Create the mesh using the projected data so all measurements are in m
# The function inla.mesh.2d() is recommended and will create
# Constrained Refined Delaunay Triangulation (CRDT) == mesh
# The boundary will have a variance 2x the variance of the domain
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
###############################################################################

# use a mesh based on a non-convex hull
# to avoid adding many small triangles outside the domain of interest
# (more triangles = larger computation times)
# inla.nonconvex.hull(points, convex, concave, resolution, eps)
# points = 2D point coordinates (2-column matrix)
# convex = desired extension radius
#   Determines the smallest allowed convex curvature radius.
#   Negative values are fractions of the approximate initial set diameter.
# concave = minimal concave curvature radius
#   Default is concave = convex.
# resolution = internal computation resolution
#   A warning will be issued when this needs to be increased for
#   higher accuracy, with the required resolution stated.
# eps = The polygonal curve simplification tolerance used for simplifying
#   the resulting boundary curve.
prdomain <- inla.nonconvex.hull(coordinates(data.all),
  convex = -0.05, resolution = c(40, 15))

# inla.mesh.2d(loc = NULL, loc.domain = NULL, offset = NULL,
#   n = NULL, boundary = NULL, interior = NULL, max.edge,
#   min.angle = NULL, cutoff = 1e-12, plot.delay = NULL)
# loc = locations used as initial triangulation nodes
# loc.domain = a polygon to determine the domain extent
# max.edge = specifies the maximum allowed triangle edge lengths
#   in the inner domain and in the outer extension.
#   A scalar or length two vector, on the SAME SCALE UNIT as
#   the coordinates.
# offset = a numeric, or length two vector.
#   If negative then it is a factor relative to the approx. data diameter.
#   If positive it is the extension distance on same scale as coords.
# n = initial number of points on the extended boundary
# interior = a list of segments to specify interior constraints,
#   each one of inla.mesh.segment class.
#   A good mesh needs to have triangles as regular as possible
#   in size and shape.
# min.angle = a scalar or length two vector,
#   to specify the minimum internal angles of the triangles on
#   the inner domain and on the outer extension.
#   Values up to 21 guarantee the convergence of the algorithm.
# cutoff = minimum allowed distance between points.
#   Points at a closer distance than the supplied value are replaced
#   by a single vertex to avoid small triangles.
#   Critical when we have some very close points, either for
#   point locations or on the domain boundary.
mesh <- inla.mesh.2d(loc = coordinates(data.all),
  offset = c(-0.06, -0.03),
  cutoff = 100,
  max.edge = c(100, 2000),
  boundary = prdomain)

spde <- inla.spde2.matern(mesh)

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
  X_xp <- matrix(1, ncol = 1, nrow = spde$n)

  data = list(
    n_i = NROW(data.all@data)
    n_x = spde$n,
    n_t = length(desired.years),
    n_p = ncol(X_xp),

    x_s = spde$idx$loc - 1,

    c_i = as.vector(data.all@data$WTCPUE),
    s_i = mapval[match(station, mapval$s.unique), "map"],
    t_i = factor(data.all@data$YEAR, levels = desired.years,
      labels = 1:length(desired.years)) - 1,

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
    Epsilon_input = matrix(rnorm(spde$n * data$n_t),
      nrow = spde$n, ncol = data$n_t),
    Omega_input = rnorm(spde$n)
    )

  obj <- MakeADFun(data = data,
    parameters = parameters,
    random = c("Epsilon_input", "Omega_input"),
    hessian = FALSE,
    DLL = my.tmb)

newtonOption(obj, smartsearch = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    lower = c(rep(-20, 2), rep(-10, 3), -0.999), #lower par bounds
    upper = c(rep(20, 2), rep(10, 3), 0.999),    #upper par bounds
    control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))
  opt[["final_gradient"]] <- obj$gr(opt$par)

  # Obtain standard errors
  Report <- obj$report()
  report <- try(sdreport(obj))

###############################################################################
#### Create summaries
###############################################################################

    Report <- obj$report()
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
stock <- rpart(Report$Omega ~ longitude)
stock.orig <- stock
splits <- rep(NA, dim(stock$cptable)[1] - 1)

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
    "Report" = Report, "mesh" = mesh,
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
