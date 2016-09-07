#' Calculate management units from a \code{information} read in using
#' \code{read_results}.
#'
#' @param information
#' @param file A character value giving the filename to save the resulting
#' plot of the groups within the population. The default is \code{NULL},
#' which leads to no plot being saved.
#' @param choose Should the algorithm with maximum management units be ran?
#' The default is \code{FALSE}.
#'
#' @examples
#' information <- read_results(data = sim_data, report = Report)
#' manunit <- calc_manunit(information = information)
#'
calc_manunit <- function(information, file = NULL, choose = FALSE) {

  # Pull out relevant portions of information
  info <- information$info

  # Find the clusters if you specify the true number of clusters that there
  # should be based on the length of alpha.
  cluster.true <- SPODT::spodt(omega ~ 1, data = info,
    min.parent = 4, min.child = 2,
    level.max = ifelse(NCOL(information$alpha) == 1, 2, NCOL(information$alpha)))
  lines.true <- SPODT::spodtSpatialLines(cluster.true, data = info)

  if (choose) {
    # Let the algorithm choose the clusters without providing a maximum
    # rtwo.min can be between 0 and 1, default is 0.001
    cluster.choose <- SPODT::spodt(omega ~ 1, data = info,
      level.max = NROW(information$mesh$loc), graft = 0,
      min.parent = 4, min.child = 2, rtwo.min = 0.05)
    # length(unique(cluster.choose@partition))
    lines.choose <- SPODT::spodtSpatialLines(cluster.choose, data = info)
  }

  grouppoly <- SpatialPolygons(
    tapply(1:length(information$group), information$group,
    function(x, data = information$loc_true) {
      X <- data[x, ]
      X.chull <- chull(X)
      X.chull <- c(X.chull, X.chull[1])
      Polygons(list(Polygon(X[X.chull, ])), parent.frame()$i[])
  }))

  if (!is.null(file)) {
    png(file, units = "in", res = 300, width = 8, height = 5)
    plot(info, cex = abs(min(info$omega)) + info$omega,
      col = ifelse(info$omega >= 0, "green", "red"), pch = 1)
    lines(lines.true, lty = 2, lwd = 2, xpd = TRUE)
    axis(side = 1, pos = 0)
    axis(side = 2, at = seq(0, 1500, by = 500))
    lines(grouppoly)
    dev.off()
  }
  # SPODT::test.spodt()
  # formula, data, R2.obs, rdist, par.rdist, nb.sim, weight = FALSE,
  # graft = 0, level.max = 5, min.parent = 10, min.child = 5,
  # rtwo.min = 0.001)
  # test <- SPODT::test.spodt(info@data$omega ~ 1, data = info,
  #   # R2.obs = 0.5,
  #   R2.obs = cluster.true@R2,
  #   rdist = "rnorm",
  #   par.rdist = c(NROW(info@data), mean(sd(info@data$omega)), sd(info@data$omega)),
  #   nb.sim = 10,
  #   level.max = NROW(information$mesh$loc),
  #   min.parent = 1, min.child = 1, weight = TRUE)

###############################################################################
# Calculate statistics based on the identified clustering
###############################################################################
  # todo: find the true group for each location on the mesh
  diffs <- data.frame(
    "loc" = info@coords,
    "true" = info@data$true,
    "spatial" = cluster.true@partition
    )
  diffs$spatial <- as.numeric(
    factor(diffs$spatial, labels = 1:length(unique(diffs$spatial))))

  diffs$diffs <- diffs$true - diffs$spatial
  match.table <- table(diffs$diffs)

  bad <- tapply(diffs$spatial, diffs$true, function(x) {
    tab <- table(x)
    if (dim(tab) > 1) {
      stat <- sum(tab[-which.max(tab)])
    } else { stat <- 0 }
    return(stat)
  })

  return(list("cluster.true" = cluster.true, "match" = match.table,
    "wronggroup" = bad, "grouppoly" = grouppoly))
}
