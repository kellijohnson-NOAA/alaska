#' Calculate management units from a \code{information} read in using
#' \code{read_results}.
#'
#' @param information
#' @param file A character value giving the filename to save the resulting
#' plot of the groups within the population. The default is \code{NULL},
#' which leads to no plot being saved.
#'
#' @examples
#' information <- read_results(data = sim_data, report = Report)
#' manunit <- calc_manunit(information = information)
#'
calc_manunit <- function(information, file = NULL) {

  # Pull out relevant portions of information
  info <- information$info
  # Find the full polygon
  projection <- proj4string(information$info)
  bounds <- calc_meshbound(information$mesh, projection = projection)

  if (is.null(information$lines_grouptrue)) {
    pol.true <- bounds$poly
  } else {
    pol.true <- calc_polys(bounds$poly, information$lines_grouptrue)
  }
  # Add a small amount of space to the polygon and then determine which
  # true group each point is in.
  info@data$true <- sp::over(info,
    rgeos::gBuffer(pol.true, width = 1, byid = TRUE))

    # Let the algorithm choose the clusters without providing a maximum
    # rtwo.min can be between 0 and 1, default is 0.001
    cluster.choose <- SPODT::spodt(omega ~ 1, data = info,
      level.max = NROW(information$mesh$loc), graft = 0,
      min.parent = 4, min.child = 2,
      # base the minimum variance of the partition on the estimated
      # marginal variance of Omega
      rtwo.min = information$SigmaO)
    # Create the lines for each management unit based on the chosen clusters
    lines.choose <- SPODT::spodtSpatialLines(cluster.choose, data = info)
    # split the full polygon by the estimated polygons
    if (length(unique(cluster.choose@partition)) == 1) {
      pol.choose <- pol.true
    } else {pol.choose <- calc_polys(bounds$poly, lines.choose)}

  if (!is.null(file)) {
    png(file, units = "in", res = 300, width = 8, height = 5)
    plot(info, cex = abs(min(info$omega)) + info$omega,
      col = ifelse(info$omega >= 0, "green", "red"), pch = 1,
      ylab = "northing")
    for (ii in 1:length(pol.choose)) {
      lines(pol.choose[ii], lty = ii + 1, lwd = 2, xpd = TRUE)
    }
    axis(1, pos = 0)
    text(x = mean(par("usr")[1:2]), y = -400, label = "easting")
    axis(2, at = seq(0, round(extent(info)@ymax / 5, -2) * 5, length.out = 4))
    dev.off()
  }

###############################################################################
# Calculate statistics based on the identified clustering
###############################################################################
  # Create a data frame that has the true groups and the estimated groups
  diffs <- data.frame(info@coords,
    "true" = info@data$true,
    "est" = cluster.choose@partition)

  match.table <- tapply(cluster.choose@partition, info@data$true, table)
  if (is.list(match.table)) {
    match.table <- do.call("rbind", match.table)

  # Make match.table a data frame and correctly set the names, which
  # are lost when using a matrix
  match.table <- data.frame(match.table)
  colnames(match.table) <- gsub("X", "", colnames(match.table))

  correctlabel <- rep(NA, NROW(cluster.choose@adj))
  # 1. Assign that number w/ most matches to pol.choose
  # 2. Remove that pol.true and pol.choose from the options
  # 3. Repeat the process until there are no more rows options
  # 4. Find which numbers were not assigned b/c they do not correspond to
  # 5. a true group and give them a number
  if (NROW(match.table) == 1) {
    best <- which(match.table == max(match.table), arr.ind = TRUE)
    correctlabel[best[1]] <- rownames(match.table)[best[1]]
  }
  while (NROW(match.table) > 1) {
    best <- which(match.table == max(match.table), arr.ind = TRUE)
    correctlabel[best[1]] <- rownames(match.table)[best[1]]
    if (NROW(match.table) == 2 & NCOL(match.table) == 2) {
      best <- which(match.table == match.table[-best[1], -best[2]],
        arr.ind = TRUE)
      correctlabel[best[1]] <- rownames(match.table)[best[1]]
    }
    match.table <- match.table[-best[1], -best[2]]
    rm(best)
  }
  while (any(is.na(correctlabel))) {
    correctlabel[which(is.na(correctlabel))[1]] <-
      seq_along(correctlabel)[!seq_along(correctlabel) %in% correctlabel][1]
  }
  info@data$est <- factor(cluster.choose@partition, labels = correctlabel)
  } else {
    info@data$est <- 1
  } # End if (is.list(match.table))

###############################################################################
# Calculate areas of overlap between groups
###############################################################################
  areas <- matrix(NA, nrow = length(pol.true), ncol = length(pol.choose),
    dimnames = list(paste0("true.", seq_along(pol.true)),
    paste0("choose.", seq_along(pol.choose))))
  for (x in seq_along(pol.true)) {
    for (y in seq_along(pol.choose)) {
      temp <- rgeos::gIntersection(pol.choose[y], pol.true[x])
      if (is.null(temp)) {
        areas[x, y] <- 0
      } else {areas[x, y] <- rgeos::gArea(temp) / rgeos::gArea(pol.true[x])}
      rm(temp)
    }
  }

  return(list("cluster.choose" = cluster.choose, "areas" = areas,
    "info" = info))
}
