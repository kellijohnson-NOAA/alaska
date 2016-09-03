###############################################################################
###############################################################################
## Purpose:  Stock analysis of cod in Alaska
##           Import shapefiles of Alaska
## Author:   Kelli Faye Johnson
## Contact:  kellifayejohnson@gmail.com
## Date:     2016-09-03
## Comments:
###############################################################################
###############################################################################

###############################################################################
#### Import shapefile
###############################################################################
efhCRS <- CRS("+proj=aea +lat_1=65.0 +lat_2=55.0
               +lat_0=50.0 +lon_0=-154.0
               +x_0=0 +y_0=0 +units=m")
maps.efh <- readShapePoly(file.path(dir.data,
                                    "efh_shapefile_2005", "EFH_2005"),
                          proj4string = efhCRS)
maps.efh <- spTransform(maps.efh, akCRS)

names.nt <- c("Near_Shore_Bristol_No_Trawl", "Northern_BS_ResearchArea",
              "Nunivak_Kusko", "Prib_Hab_Cons_Area", "St_Lawerence", "St_Matts",
              "Modified_Gear_Trawl_Zone", "SE_No_Trawl", "Cook_inlet",
              "GOA_Type_1", "AI_HCA", "BS_HCA", "GOA_Slope_HCA",
              "Red_King_Crab_Closure_Area")
names.nf <- c("Bowers_Ridge", "AK_Seamount_HPA")

for(q in seq_along(names.nt)){
  assign(paste0("names.nt", "_", q),
    readShapePoly(file.path(dir.data, "alaska_SSLShapefiles", names.nt[q]),
    proj4string = efhCRS))
}
for(q in seq_along(names.nf)){
  assign(paste0("names.nf", "_", q),
    readShapePoly(file.path(dir.data,
    "alaska_SSLShapefiles", names.nf[q]), proj4string = efhCRS))
}
maps.eez <- readShapeSpatial(file.path(dir.data, "USMaritimeLimitsAndBoundaries",
  "USMaritimeAlaskaEEZ"), proj4string = efhCRS)
maps.ak <- readShapeSpatial(file.path(dir.data, "USMaritimeLimitsAndBoundaries",
  "alaska_coast"), proj4string = efhCRS)
