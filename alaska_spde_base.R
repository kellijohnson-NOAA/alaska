###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2015-01-05
## Input:      All data was downloaded from
##             http://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm
##             and placed in dir.data
## Comments:
###############################################################################
###############################################################################

###############################################################################
#### Set initial inputs
###############################################################################
# Must set my.base or have loaded in R
my.pars      <- "par(oma = rep(0, 4), mar = c(2, 2.7, 0, 0.5))"
my.tmb       <- "gompertz_kfj"

my.filetype <- "png"
my.width <- c(3.34646, 6.69291) # 85 mm, 170mm
my.height <- my.width[1]
my.height.map <- 4.5
my.resolution <- 500

getweb <- FALSE #If TRUE pull r4kfj off of github
run.all <- TRUE #logical; save or load(saved) workspace

desired.areas <- c("ai", "goa")
desired.spp <- "Pacific cod"
desired.years <- 1990:2015 #2015 = last year of data in the AI

###############################################################################
#### Load packages
###############################################################################
    dir.results <- file.path(my.base, "results")
    dir.data    <- file.path(my.base, "data")
    dir.create(dir.results, showWarnings = FALSE)
    dir.create(dir.data, showWarnings = FALSE)

  # Install personal R package: "r4kfj"
  if(getweb == TRUE) {
    devtools::install_github("r4kfj", "kellijohnson")
  }

  # Load .R files specific for the alaska analysis, located in "lib" folder
  sapply(dir(file.path(my.base, "lib"), full.names = TRUE), source)
  library("r4kfj")
  load_packages(c("igraph", "INLA", "maps", "maptools", "mapproj", "Matrix",
                  "plyr", "raster", "rpart", "reshape", "reshape2", "rgdal",
                  "rgeos", "sos", "sp", "spdep", "stats4", "TMB", "xtable"))

  define_projections()
