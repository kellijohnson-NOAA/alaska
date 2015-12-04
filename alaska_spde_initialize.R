###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2015-01-05
## Comments:   Polygons for the seven areas in Alaska are based on a shape file
##             provided by Annie Greig of NOAA (email: angie.greig@noaa.gov).
##             Shape file were in arcgis format.
##             The function create_areas assumes the shapefile is pocentered 
##             or projected (i.e. not in latlon format)
###############################################################################
###############################################################################

###############################################################################
#### Get data
###############################################################################
if(run.all == TRUE){
  # create a vector of file names containing needed data
  # Read in all of the raw data in the dir.data
    race.num <- create_race_num(speciesnames = desired.spp)
    
    data.files <- unlist(sapply(desired.areas, function(x) {
                                dir(file.path(dir.data, "all.data"), 
                                    full.names = TRUE, pattern = x)
                                }
                                )
                        )
  
  # Remove entries without station numbers & Japanese trawl data
  # Read in all files in the dir.data folder that contain 
  # the areas listed in desired.areas
    data.all <- subset(do.call("rbind",
                               lapply(data.files, read.csv, 
                                      na.strings = -9999)),
                       STATION != "" & VESSEL < 500)

    data.all$station <- with(data.all, paste(STATION, STRATUM, sep = "_"))
    data.all$id <- with(data.all, paste(station, SID, sep = "_"))
  # Add coordinate data to the data.all, using the akCRS projection
    coordinates(data.all) <- with(data.all, cbind("LONGITUDE", "LATITUDE"))
    proj4string(data.all) <- llCRS
    data.all <- spTransform(data.all, akCRS)

  # GIS Analysis
  # shape.manage <- create_areas(shape = shape.afsc)

  # Polygons for areas, that cuts off the eastern portion of the gulf of AK
  # because of the management area.
  alaska_areas_management <- list(
      "GOA_W" = matrix(c(c(-140, -170, -170, -140, -140) %% 360,
                         c(45.0, 45.0, 65.0, 65.0, 45.0)), ncol = 2),
      "AI" = matrix(c(c(165, 165, -170, -170, 165) %% 360,
                          c(45.0, 65.0, 65.0, 45.0, 45.0)), ncol = 2)
      )
  spatialareas <- SpatialPolygons(lapply(alaska_areas_management, function(x) {
        Polygons(list(Polygon(x)), 
                 ID = names(alaska_areas_management)[substitute(x)[[3]]])
        }))
    proj4string(spatialareas) <- llCRS
    spatialareas <- spTransform(spatialareas, akCRS)

  data.all$inside <- over(data.all, spatialareas)

  # Subset the data for the years and values inside the study area
  data.all <- subset(data.all, !is.na(inside) & YEAR %in% desired.years)

  # Subset data for those entries where the desired spp was not found
  # This will be useful for rerunning the model with a small value added
  # to the tows where the species was not caught to test the sensitivity
  # to the assumption that I do not need to model the zero tow data
  data.zero <- subset(data.all, !SID %in% race.num$RACE, 
                      select = c("STATION", "STRATUM", "YEAR",
                                "DATETIME", "WTCPUE", 
                                "SID", "station", "id", "inside"))
  temp <- with(data.zero@data, paste(station, YEAR))
  data.zero.keep <- sapply(split(seq_along(temp), temp), "[", 1)
  temp <- list()
  for(spp in seq_along(desired.spp)) {
    temp[[spp]] <- data.zero[data.zero.keep, ]
    temp[[spp]]$SID <- race.num$RACE[spp]
  }
  data.zero <- do.call("rbind", temp)
  
  # Subset data 
  data.spp <- subset(data.all, SID %in% race.num$RACE, 
                     select = c("STATION", "STRATUM", "YEAR",
                                "DATETIME", "WTCPUE", 
                                "SID", "station", "id", "inside"))

  save.image(file.path(dir.results, "alaska_spde_initialize.RData"))
} else {
  load(file.path(dir.results, "alaska_spde_initialize.RData"))
}
