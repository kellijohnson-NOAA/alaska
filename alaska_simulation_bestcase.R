###############################################################################
###############################################################################
## Purpose:    Simulate data for a spatial gompertz model
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-08-30
## Comments:
###############################################################################
###############################################################################
if(basename(getwd()) != "alaska") stop("set your working directory",
  " to the cloned repo")
library(INLA)
library(TMB)
dir.results <- file.path(getwd(), "results")
dir.create(dir.results)
sapply(dir("lib", full.names = TRUE), source)
tmb <- gsub("\\.cpp", "", dir("data", pattern = "cpp"))[1]
if(!file.exists(file.path("data", paste0(tmb, ".o")))) {
  TMB::compile(file.path("data", paste0(tmb, ".cpp")))
}
dyn.load(dynlib(file.path("data", tmb)))

###############################################################################
# Simulate data
###############################################################################
data_i <- Sim_Gompertz_Fn(
  n_years = 50, #50,
  n_stations = 500, #500,
  phi = 0.0,
  SpatialScale = 0.01,
  SD_O = 0.2,
  SD_E = 0.1,
  SD_obs = 0.01,
  rho = 0.5,
  logMeanDens = 1.0,
  Loc = NULL,
  projection = NULL,
  seed = 10)

###############################################################################
#### Build inputs
###############################################################################
mesh_i <- INLA::inla.mesh.create(data_i[, c("Longitude", "Latitude")])
obj <- calc_adfun(data = data_i, mesh = mesh_i, 
  tmb = tmb,
  variable = "zeroinflatedlnorm")

# Run optimizer
optimizer <- nlminb(obj$par, objective = obj$fn,
  gradient = obj$gr,
  lower = c(rep(-200, 5), -0.999, rep(-200, 2)),
  upper = c(rep( 200, 5), 0.999, rep(200, 2)),
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1))

save(data_i, obj, file = file.path("results", "simulation_bestcase.RData"))
sink(file.path("results", "simulation_bestcase.txt"))
cat("Results from the best case simulation\n")
cat("exp(theta_z)\n")
cat(exp(obj$report()$theta_z))
cat("\n\n")
unlist(obj$report()[c("SigmaO","SigmaE", "Range", "rho", "alpha", "phi")])
sink()
plot(data_i$Omega[data_i$Year == 1], obj$report()$Omega_x[1:length(unique(data_i$Site))]); abline(0, 1)
test <- matrix(data_i$Epsilon, nrow = length(unique(data_i$Site)), byrow = TRUE)
test2 <- matrix(obj$report()$Epsilon_xt, ncol = length(unique(data_i$Year)))
plot(test, test2[1:NROW(test), ]); abline(0, 1)
plot(obj$report()$c_i, exp(obj$report()$log_chat_i)); abline(0, 1)

rm(data_i, obj, mesh_i, test, test2, tmb)
# End of file
