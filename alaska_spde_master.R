###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2014-08-17
## Comments:   Polygons for the areas in Alaska are based on a shape file
##             provided by Annie Greig of NOAA (email: angie.greig@noaa.gov).
##             Shape file were in arcgis format.
###############################################################################
###############################################################################

###############################################################################
#### Set initial inputs
###############################################################################
my.base <- file.path("c:", "alaska")
  if (!file.exists(my.base)) {
    my.base <- file.path("d:", "alaska")
    if(!file.exists(my.base)) stop(paste(my.base, "does not exist"))
  }
setwd(my.base)

getweb <- FALSE #If TRUE pull r4kfj off of github

###############################################################################
#### Run analysis
###############################################################################
source("alaska_spde_initialize.R")
source("alaska_spde_data.R")
source("alaska_spde_mesh.R")
source("alaska_shapefiles.R")

# Initial figure
source("alaska_alaskamap.R")
source("alaska_alaskadata.R")
source("alaska_fig_data.R")
source("alaska_fig_mesh.R")

source("alaska_cleanup.R")

# Simulations
source("alaska_simulation_SpatialScale.R")

# Analysis
source("alaska_spde_analysis.R")
source("alaska_spde_results.R")

# Final clean up
source("alaska_cleanup02.R")
