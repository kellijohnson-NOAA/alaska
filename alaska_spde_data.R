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
# Get data
###############################################################################
# Eventually data will be subset for these columns
keepcolumns <- c("STATION", "STRATUM", "YEAR", "DATETIME", "WTCPUE",
  "SID", "station", "id", "inside")

# create a vector of file names containing needed data
# Read in all of the raw data in the dir.data
race.num <- read.csv(file.path(dir.data, "EBS_RACE_Look_2012.csv"))
race.num <- data.frame(race.num[race.num$COMMON %in% desired.spp, "RACE"])
colnames(race.num) <- "RACE"
row.names(race.num) <- desired.spp

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

# Assign all other species to a weight of zero
# Remove duplicate stations such that the final data set has
# both positive and zero tows for the species of interest
data.all$WTCPUE[!data.all$SID %in% race.num$RACE] <- 0
data.all$NUMCPUE[!data.all$SID %in% race.num$RACE] <- 0
data.all$SID[!data.all$SID %in% race.num$RACE] <- race.num$RACE + 1
data.all$station <- with(data.all, paste(STATION, STRATUM, sep = "_"))
data.all$id <- with(data.all, paste(station, SID, sep = "_"))
data.all <- data.all[order(data.all$YEAR, data.all$DATETIME, data.all$id), ]
data.all <- data.all[!duplicated(with(data.all, paste(YEAR, station))), ]

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
