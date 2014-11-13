###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create spde for spatial analysis
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2014-08-17
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
  shape.afsc <- readOGR(file.path(dir.data, "stratum"), "strata_geographic")
  shape.afsc.akcrs <- spTransform(shape.afsc, akCRS)
  shape.manage <- create_areas(shape = shape.afsc.akcrs)
  data.all$inside <- over(data.all, 
                          as(eval(parse(text = my.shape)), "SpatialPolygons"))

  # Subset data 
  data.spp <- subset(data.all, !is.na(inside) & 
                     SID %in% race.num$RACE & YEAR %in% desired.years, 
                     select = c("STATION", "STRATUM", "YEAR",
                                "DATETIME", "WTCPUE", 
                                "SID", "station", "id", "inside"))

  save.image(file.path(dir.results, "alaska_spde_initialize.RData"))
} else {
  load(file.path(dir.results, "alaska_spde_initialize.RData"))
}
