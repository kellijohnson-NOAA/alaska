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

create_areas <- function(shape = shape.afsc.akcrs,
                         listofareas = alaska_areas,
                         subset = NA){
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

    listofpolygons <- lapply(listofareas, function(x) {
            Polygons(list(Polygon(x)), ID = names(listofareas)[substitute(x)[[3]]])
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
    if (any(!is.na(subset))) {
        return(shape.full[which(names(shape.full) %in% subset)])
    } else {
        return(shape.full)
    }
}

# endoffile