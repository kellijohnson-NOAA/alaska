#' Plot a general map of the study area
#'
#' Description 
#' @author Kelli Faye Johnson
#' @param map.shape A SpatialPointsPolygon to be ploted as the base map
#' @param axis Logical to plot axis labels or not
#' @param labels.major Logical to add major labels
#' @param labels.canyons Logical to add canyon labels
#' @param labels.current Logical to add currents and
#' their associated arrows for the current directions
#' @param text.col Text colour
#' @export

    map_alaska <- function(map.shape = map.outline, axis = TRUE, 
                           text.col = "darkred", font = 2, 
                           labels.major = TRUE, labels.minor = FALSE,
                           labels.canyons = FALSE, labels.currents = FALSE,
                           ...) {
        plot(map.shape, ...)
        if (axis) {
          r4kfj::llgridlines(map.shape, recenter = TRUE)
        }
        places.major <- data.frame("name" = c("Alaska", "Bering Sea", 
                                            "Gulf of\nAlaska", "Aleutian Islands"),
                                   "LATITUDE" = c(62, 59, 58, 50),
                                   "LONGITUDE" = c(203, 190, 215, 190),
                                   "cex" = rep(1, 4),
                                   "pos" = rep(1, 4),
                                   "offset" = rep(0, 4)
                                   )
        places.minor <- data.frame("name" = c("Dixon \nEntrance"),
                                   "LATITUDE" = c(54),
                                   "LONGITUDE" = c(227),
                                   "cex" = rep(0.6, 1),
                                   "pos" = rep(1, 1),
                                   "offset" = rep(0, 1)
                                   ) 
        places.canyons <- data.frame("name" = c("Unimak\nPass", "Samalga\nPass", 
                                              "Amchitka\nPass", "Buldir\nStrait"),
                                   "LATITUDE" = c(54.4, 52.5, 51.2, 52),
                                   "LONGITUDE" = c(197, 190, 180, 173),
                                   "cex" = rep(0.6, 4),
                                   "pos" = rep(1, 4),
                                   "offset" = rep(0.5, 4)
                                   )
        places.currents <- data.frame("name" = c("Alaska\nCoastal\nCurrent", 
                                                 "Alaskan Stream", 
                                                 "North Slope Current"),
                                      "LATITUDE" = c(55, 53.5, 55),
                                      "LONGITUDE" = c(226, 211, 190),
                                      "cex" = rep(0.6, 3),
                                      "pos" = rep(2, 3),
                                      "offset" = rep(0.5, 3)
                                      )
        places.arrows <- data.frame("LONGITUDE" = c(227.2, 207.3, 199.0, 196, 192.6, 212, 199.8, 180.9, 180, 175.7, 175.6,
                                                    208, 201.3, 194.3, 192.6, 190.3, 200, 181, 180.2, 176, 171, 191.5),
                                    "LATITUDE" = c(54.5, 58.7, 55.2, 54.0, 53.0, 57.2, 53, 51, 51, 51.5, 53,
                                                   58.7, 55.5, 54.6, 53.0, 53.4, 53, 51, 52, 51.5, 53, 54),
                                    "curve" = rep(c(-0.47, 0, 0.5, 0, .76, 0.2, 0.1, 0.6, 0.1, 0.3, -0.2), 2))
       fxn.env <- new.env()
       places.akCRS <- sapply(paste("places", c("major", "minor", "canyons", "currents", "arrows"), sep = "."), function(x) {
          data <- eval(as.name(x))
          coordinates(data) <- data[, c("LONGITUDE", "LATITUDE")]
          proj4string(data) <- llCRS
          data <- spTransform(data, CRS(proj4string(map.shape)))
          assign(x, data, envir = fxn.env)
        })
        if(labels.major) {
            text(fxn.env$places.major, labels = places.major$name, 
                 col = text.col, cex = places.major$cex,
                 font = font, offset = places.major$offset, 
                 pos = places.major$pos)
        }
        if(labels.minor) {
            text(fxn.env$places.minor, labels = places.minor$name, 
                 col = text.col, cex = places.minor$cex,
                 font = font, offset = places.minor$offset, 
                 pos = places.minor$pos)
        }        
        if(labels.canyons) {
            text(fxn.env$places.canyons, labels = places.canyons$name,
                 col = text.col, cex = places.canyons$cex,
                 font = font, offset = places.canyons$offset, 
                 pos = places.canyons$pos)
        }
        if(labels.currents) {
            text(fxn.env$places.currents, labels = places.currents$name,
                 col = text.col, cex = places.currents$cex,
                 font = font, offset = places.currents$offset, 
                 pos = places.currents$pos)
            places.arrows.num <- dim(fxn.env$places.arrows)[1] / 2
            places.arrows.plot <- cbind(coordinates(fxn.env$places.arrows[1:places.arrows.num, ]),
                                        coordinates(fxn.env$places.arrows[(places.arrows.num + 1):(places.arrows.num * 2), ]),
                                        fxn.env$places.arrows[1:places.arrows.num, "curve"]@data)
            apply(places.arrows.plot, 1, function(x) {
                igraph:::igraph.Arrows(x[1], x[2], x[3], x[4],
                              h.lwd = 1.5, sh.lwd = 1.5,
                              curve = x[5], width = 1, size = 0.4)
                })
            }
        }
