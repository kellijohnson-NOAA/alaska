###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-04-13
## Comments:   Alaska polygons are based on a shape file provided by
##             Annie Greig of NOAA (email: angie.greig@noaa.gov).
##             Shape file were in arcgis format.
##             http://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm
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

  if (getweb) {
    # Install personal R package: "r4kfj"
    devtools::install_github("kellijohnson/r4kfj@master")
    # Install the stable version of INLA
    install.packages("INLA", repos = "http://www.math.ntnu.no/inla/R/stable")
    # Install a needed package for rgeos
    install.packages("gpclib", repos = "http://cran.fhcrc.org/", type = "source")
  }
  # Install TMB
  if (file.exists(file.path("c:", "adcomp"))) {
    old_wd <- getwd()
    setwd(file.path("c:", "adcomp"))
    if (getweb) {
      system("git fetch")
      system("git rebase origin/master")
    }
    remove.packages("TMB")
    source("install_windows.R")
    setwd(old_wd)
    rm(old_wd)
  } else {
    stop("TMB is not installed on this computer.")
  }
  # Load .R files specific for the alaska analysis, located in "lib" folder
  ignore <- sapply(dir(file.path(my.base, "lib"), full.names = TRUE), source)
  library("r4kfj")
  library(INLA)
  load_packages(c("igraph", "maps", "maptools", "mapproj", "Matrix",
                  "plyr", "raster", "rpart", "reshape", "reshape2", "rgdal",
                  "rgeos", "sos", "sp", "splancs",
                  "spdep", "stats4", "TMB", "xtable"))

  define_projections()

###############################################################################
#### Get data
###############################################################################
  # Eventually data will be subset for these columns
  keepcolumns <- c("STATION", "STRATUM", "YEAR", "DATETIME", "WTCPUE",
    "SID", "station", "id", "inside")

  # create a vector of file names containing needed data
  # Read in all of the raw data in the dir.data
  race.num <- create_race_num(speciesnames = desired.spp)

  data.files <- dir(file.path(dir.data, "all.data"),
    pattern = paste0(desired.areas, collapse = "|"), full.names = TRUE)

  # Read in all files in the dir.data folder that contain
  # the areas listed in desired.areas
  # Remove entries without station numbers & Japanese trawl data
  # Remove early years
  data.all <- do.call("rbind",
    lapply(data.files, read.csv, na.strings = -9999))
  data.all <- data.all[data.all$YEAR %in% desired.years, ]
  data.all <- data.all[data.all$VESSEL < 500, ]
  data.all <- data.all[data.all$STATION != "", ]
  data.all$station <- with(data.all, paste(STATION, STRATUM, sep = "_"))
  data.all$id <- with(data.all, paste(station, SID, sep = "_"))

  # Add coordinate data to the data.all, using the akCRS projection
  coordinates(data.all) <- with(data.all, cbind("LONGITUDE", "LATITUDE"))
  proj4string(data.all) <- llCRS
  data.all <- spTransform(data.all, akCRS)

  # Polygons for areas, that cuts off the eastern portion of the gulf of AK
  # because of the management area.
  alaska_areas_management <- list(
    "GOA_W" = matrix(c(c(-140, -170, -170, -140, -140) %% 360,
                       c(45.0, 45.0, 65.0, 65.0, 45.0)), ncol = 2),
    "AI" = matrix(c(c(165, 165, -170, -170, 165) %% 360,
                    c(45.0, 65.0, 65.0, 45.0, 45.0)), ncol = 2)
      )
  spatialareas <- SpatialPolygons(sapply(seq_along(alaska_areas_management),
    function(x, y = alaska_areas_management) {
    Polygons(list(Polygon(y[x])), ID = names(y)[x])
    }))
    proj4string(spatialareas) <- llCRS
    spatialareas <- spTransform(spatialareas, akCRS)

  # Subset the data for the years and values inside the study area
  data.all$inside <- over(data.all, spatialareas)
  data.all <- data.all[!is.na(data.all$inside), ]

  # Subset data for those entries where the desired spp was not found
  # This will be useful for rerunning the model with a small value added
  # to the tows where the species was not caught to test the sensitivity
  # to the assumption that I do not need to model the zero tow data
  data.zero <- list()
  for (spp in seq_along(desired.spp)) {
    data.zero[[spp]] <- data.all[!data.all$SID %in% race.num$RACE[spp],
      keepcolumns]
    data.zero[[spp]] <- data.zero[[spp]][
      !duplicated(data.zero[[spp]]@data[, c("STATION", "STRATUM", "YEAR")]), ]
    data.zero[[spp]]@data$WTCPUE <- 0
    data.zero[[spp]]@data$SID <- race.num$RACE[spp]
  }
  data.zero <- do.call("rbind", data.zero)

  # Subset data
  data.spp <- data.all[data.all$SID %in% race.num$RACE, keepcolumns]
