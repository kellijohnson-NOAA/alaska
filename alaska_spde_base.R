###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2014-08-17
## Comments:   
###############################################################################
###############################################################################

###############################################################################
#### Set initial inputs
###############################################################################
my.data.name    <- "data.spp"
my.pars         <- "par(oma = rep(0, 4), mar = c(2, 2.7, 0, 0.5))"
my.shape        <- "shape.manage"
my.tmb          <- "gompertz_kfj"

my.filetype <- "png"
my.width <- c(3.34646, 6.69291) # 85 mm, 170mm
my.height <- my.width[1]
my.height.map <- 4.5
my.resolution <- 500

getweb <- FALSE #If TRUE pull r4kfj off of github
run.all <- TRUE #logical; save or load(saved) workspace

desired.areas <- c("ai", "goa")
desired.spp <- "Pacific cod"
desired.years <- 1990:2013 #2013 = last year of data, 2013 survey in BS & GOA


###############################################################################
#### Load packages
###############################################################################
  setwd(my.base)
    dir.results <- file.path(my.base, "results")
    dir.data    <- file.path(my.base, "data")

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