#'
#'
calc_meshmap <- function(mesh) {
  #A.locator is the location of each data point on the mesh
  #Only some of the mesh nodes are actually filled
  #Filled nodes == A.locator.unique
  #station is used to map 1:numberofnodes to 0:numberofnodes-1
  #station_map creates a vector to map the 0:numberofusednodes-1
  #to 0:numberofnodes-1 because the cpp code only predicts for used nodes
  A.locator <- mesh$idx$loc
    A.locator.unique <- unique(A.locator)[order(unique(A.locator))]
    station <- A.locator - 1
    station_unique  <- A.locator.unique - 1
    mapval <- data.frame(s.unique = station_unique,
                         map = seq(0, (length(A.locator.unique) - 1)))
    final <- mapval[match(station, mapval$s.unique), "map"]
    return(final)
}
