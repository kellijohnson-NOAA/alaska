#' Calculate management units from a \code{readin} read in using
#' \code{read_results}.
#'
#' @param readin A list of information read in from the disk using
#' \code\link{read_results} that summarizes the data used to simulate the
#' population and the model results.
#' @param file A character value giving the filename to save the resulting
#' plot of the groups within the population. The default is \code{NULL},
#' which leads to no plot being saved. A full path can be specified, or else
#' the file will be saved in the current directory.
#'
#' @examples
#' temp <- read_results(data = sim_data, report = Report)
#' manunit <- calc_manunit(readin = temp)
#' rm(temp, manunit)
#'
calc_manunit <- function(readin, file = NULL) {
  # Pull out relevant portions of readin
  info <- readin$info
  # Find the full polygon
  projection <- proj4string(readin$info)
  bounds <- calc_meshbound(readin$mesh, projection = projection)

  if (is.null(readin$lines_grouptrue)) {
    pol.true <- bounds$poly
  } else {
    pol.true <- calc_polys(bounds$poly, readin$lines_grouptrue)
  }
  # Add a small amount of space to the polygon and then determine which
  # true group each point is in.
  info@data$true <- sp::over(info,
    rgeos::gBuffer(pol.true, width = 1, byid = TRUE))

    # Let the algorithm choose the clusters without providing a maximum
    # rtwo.min can be between 0 and 1, default is 0.001
    cluster.choose <- SPODT::spodt(omega ~ 1, data = info,
      level.max = NROW(readin$mesh$loc),
      graft = readin$SigmaO,
      min.parent = 4, min.child = 2,
      # base the minimum variance of the partition on the estimated
      # marginal variance of Omega
      rtwo.min = 0.001)
    # Create the lines for each management unit based on the chosen clusters
    lines.choose <- SPODT::spodtSpatialLines(cluster.choose, data = info)
    # split the full polygon by the estimated polygons
    if (length(unique(cluster.choose@partition)) == 1) {
      pol.choose <- bounds$poly
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
  match.table <- tapply(cluster.choose@partition, info@data$true, table)
  if (is.list(match.table)) {
    # Find the number of matches to each estimated partition
    match.table <- Reduce(function(...) merge(..., by = "Var1", all = TRUE),
      lapply(match.table, function(x) data.frame(x)))

    # Create a blank vector, with a single entry per each cluster.choose
    correctlabel <- rep(NA, NROW(cluster.choose@adj))
    names(correctlabel) <- match.table[, 1]
    rownames(match.table) <- match.table[, 1]

   # Remove the names and rename the columns
    match.table[, 1] <- 0
    colnames(match.table) <- 0:(NCOL(match.table) -1)

    while (NCOL(match.table) > 1) {
      col <- which.max(sapply(apply(match.table, 2, which.max),
        function(x) match.table[x, parent.frame()$i[]]))
      row <- which.max(as.numeric(as.character(match.table[, names(col)])))

      correctlabel[row] <- as.numeric(names(col))
      match.table[row, ] <- 0
      match.table <- match.table[, -col]
    }

    # Fill in missing values
    temp <- 1:length(correctlabel)
    correctlabel[which(is.na(correctlabel))] <- temp[!temp %in% correctlabel]
    rm(temp)

  info@data$est <- factor(cluster.choose@partition, labels = correctlabel)
  } else {
    info@data$est <- 1
    correctlabel <- c("1" = 1)
  } # End if (is.list(match.table))

###############################################################################
# Calculate areas of overlap between groups
###############################################################################
  areas <- matrix(NA, nrow = length(pol.true), ncol = length(pol.choose),
    dimnames = list(paste0("true.", seq_along(pol.true)),
    paste0("choose.", seq_along(pol.choose))))
  for (ix in seq_along(pol.true)) {
    for (iy in seq_along(pol.choose)) {
      temp <- rgeos::gIntersection(pol.choose[iy], pol.true[ix])
      if (is.null(temp)) {
        areas[ix, iy] <- 0
      } else {areas[ix, iy] <- rgeos::gArea(temp) / rgeos::gArea(pol.true[ix])}
      rm(temp)
    }
  }

  return(list(
    "areas" = areas,
    "cluster.choose" = cluster.choose,
    "lines.choose" = lines.choose,
    "pol.choose" = pol.choose,
    "info" = info
  ))
}

export_prm("d:/alaska/data/normal.prm", getwd(), "test")
export_geo(best[[51]]$info, getwd(), "test")
system(paste(my.SaTScan, "test.prm"))
temp <- readShapePoly(file.path(getwd(), "test.col"))
projection(temp) <- llCRS

plot(best[[int]]$info)
lines(best[[int]]$lines_grouptrue)
plot(spTransform(temp, akCRS), add = TRUE)
