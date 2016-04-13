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
##             The function create_areas assumes the shapefile is pocentered
##             or projected (i.e. not in latlon format)
###############################################################################
###############################################################################

###############################################################################
#### Set initial inputs
###############################################################################
my.base <- file.path("c:", "alaska")
  if (!file.exists(my.base)) {
    my.base <- file.path("d:", "alaska")
    if(file.exists(my.base)) stop(paste(my.base, "does not exist"))
  }
setwd(my.base)

getweb <- TRUE #If TRUE pull r4kfj off of github

###############################################################################
#### Run analysis
###############################################################################
source("alaska_spde_initialize.R")
source("alaska_spde_analysis.R")

source("alaska_spde_results.R")

