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
setwd("c:/alaska")
source("alaska_spde_base.R")
source("alaska_spde_initialize.R")
source("alaska_spde_analysis.R")

source("alaska_spde_results.R")

