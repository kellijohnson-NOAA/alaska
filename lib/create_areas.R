#' Create polygons in Alaska based on management areas
#'

#' @author Kelli Faye Johnson
#' @param shape The shapefile, must either be projected or 
#' centered on the Pacific Ocean
#' @param listofareas A list of areas to break up the Gulf of Alaska
#' and the Aleutian Islands
#' @param subset which areas (a vector of numbers) should be 
#' returned in the final SpatialPolygons item.
#' @export

create_areas <- function(shape = shape.afsc.akcrs){
    #Polygons for areas
    alaska_areas_management <- list(
      "GOA_W" = matrix(c(c(-140, -170, -170, -140, -140) %% 360,
                         c(45.0, 45.0, 65.0, 65.0, 45.0)), ncol = 2),
      "AI" = matrix(c(c(165, 165, -170, -170, 165) %% 360,
                          c(45.0, 65.0, 65.0, 45.0, 45.0)), ncol = 2)
      )
    alaska_noSE <- list(
      "GOA_W" = matrix(c(c(-140, -170, -170, -140, -140) %% 360,
                         c(45.0, 45.0, 65.0, 65.0, 45.0)), ncol = 2),
      "buldir" = matrix(c(c(174, 165, 165, 174, 175.3, 174.5, 174) %% 360,
                          c(45.0, 45.0, 65.0, 65.0, 52.5, 52.1, 45.0)), ncol = 2),
      "amchitka" = matrix(c(c(174, 174.5, 175.3, 174, 178.4, -178.68, 180, 180, 174) %% 360,
                            c(45.0, 52.1, 52.5, 65.0, 65.0, 52.82, 51.7, 45.0, 45.0)), ncol = 2),
      "samalga" = matrix(c(c(180, 180, -178.68, 178.4, -170, -170, 180) %% 360,
                           c(45.0, 51.7, 52.82, 65.0, 65.0, 45.0, 45.0)), ncol = 2)
      )
    alaska_areas_closure <- list(
      "GOA_W" = matrix(c(c(-140, -130, -130, -140, -140) %% 360,
                         c(45.0, 45.0, 65.0, 65.0, 45.0)), ncol = 2)
      )

    shape.data <- shape@data
    id.main <- shape.data$OBJECTID
    bs <- which(with(shape.data, AI_STRATA_ + GOA_STRATA) == 0)
    ai <- which(with(shape.data, AI_STRATA_ > 0))
    goa<- which(with(shape.data, GOA_STRATA > 0))
    id.main[bs] <- "BS"
    id.main[ai] <- "AI"
    id.main[goa]<- "GOA"
    
    shape.3areas <- unionSpatialPolygons(shape, id.main)
    shape.chain  <- unionSpatialPolygons(shape.3areas[c(1,3)],
                                         IDs = c(1,1))

    listofpolygons <- lapply(alaska_areas_management, function(x) {
            Polygons(list(Polygon(x)), ID = names(alaska_areas_management)[substitute(x)[[3]]])
        })

    spatialareas <- SpatialPolygons(listofpolygons)

    proj4string(spatialareas) <- llCRS
    spatialareas <- spTransform(spatialareas, CRS(proj4string(shape)))
    shape.chainareas <- gIntersection(shape.chain, spatialareas, byid = TRUE)
    shape.full <- SpatialPolygons(c(slot(shape.chainareas, "polygons"),
                                    slot(shape.3areas["BS"], "polygons")))
    proj4string(shape.full) <- CRS(proj4string(shape))
    id.names <- c(names(spatialareas), "BS")
    for(q in seq_along(id.names)) {
        slot(shape.full@polygons[[q]], "ID") <- id.names[q]
    }
     return(shape.full[which(names(shape.full) %in% 
                             names(alaska_areas_management))])
 }

# endoffile