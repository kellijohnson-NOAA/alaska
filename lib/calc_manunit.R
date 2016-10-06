#' Calculate management units from a \code{readin} read in using
#' \code{read_results}.
#'
#' @param readin A list of information read in from the disk using
#' \code\link{read_results} that summarizes the data used to simulate the
#' population and the model results.
#' @param dir The directory to save the files to. The default is \code{getwd()},
#' or your current working directory.
#' @param file A character value giving the filename to save the resulting
#' files to. No extension is necessary, as multiple files will be saved with
#' various extensions using the same name.
#' @param exe A full file path specifying the location of the SaTScan executible.
#' The default is \code{NULL}, causing \code{calc_manunit} to use
#' the \code{\link[SPODT]{spodt}}. Must use backslahes, such that the call to
#' \code{system} will work. For example: \code{"c:\\SaTScan\\SatScanBatch.exe"}.
#' @param prmtemplate A full file path specifying the location of the prm template
#' file to alter for this analysis.
#'
#' @examples
#' # Only works on local machine
#' # prm <- "d:/alaska/data/normal.prm"
#' # exe <- "c:\\SaTScan\\SaTScanBatch.exe"
#' # temp <- read_results(data = sim_data, report = Report)
#' # manunit <- calc_manunit(readin = temp, dir = getwd(),
#' #   file = "try", exe = exe, prmtemplate = prm)
#' # rm(exe, manunit, prm, temp)
#'
calc_manunit <- function(readin, dir = getwd(), file = NULL,
  exe = NULL, prmtemplate = NULL,
  longlatproj = sp::CRS("+proj=longlat +ellps=WGS84")) {
  type <- "satscan"
  if (is.null(exe)) type <- "spodt"

  #1. Pull out the relevant portions of readin
  # info is a subset of the model results for Omega,
  # that specify the point locations and the estimated values
  # that lie within the inner boundary of the mesh.
  # These were determined using \code\link{{calc_meshbound}}
  info <- readin$info
  projection <- sp::CRS(raster::projection(readin$info))

  #2. Calculate the true polygons
  bounds <- calc_meshbound(readin$mesh, projection = projection)
  if (is.null(readin$lines_grouptrue)) {
    # lines_grouptrue will be null if there is only one true subpopulation
    pol.true <- bounds$poly
  } else {
    pol.true <- calc_polys(bounds$poly, readin$lines_grouptrue)
  }
  # Add a small amount of space to the polygon and then determine which
  # true group each point is in.
  info@data$true <- sp::over(info,
    rgeos::gBuffer(pol.true, width = 1, byid = TRUE))

  #3. Estimate the management units based on spatially-explicit results
  if (type == "spodt") {
    cluster.choose <- SPODT::spodt(omega ~ 1, data = info,
      level.max = NROW(readin$mesh$loc),
      min.parent = 4, min.child = 2,
      graft = readin$SigmaO, rtwo.min = 0.001)
    # Create the lines for each management unit based on the chosen clusters
    lines.choose <- SPODT::spodtSpatialLines(cluster.choose, data = info)
    # split the full polygon by the estimated polygons
    if (length(unique(cluster.choose@partition)) == 1) {
      pol.choose <- bounds$poly
    } else {pol.choose <- calc_polys(bounds$poly, lines.choose)}
  }

  if (type == "satscan") {
    export_prm(file_in = prmtemplate, dir = dir, file = file)
    export_geo(points = info, dir = dir, file = file)
    # Double check the direction of slashes and call system
    system(paste(
      gsub("/", "\\\\", exe),
      file.path(gsub("/", "\\\\", dir), paste0(file, ".prm"), fsep = "\\")
    ), show.output.on.console = FALSE)

    pol.choose <- readShapePoly(file.path(dir, paste0(file, ".col")))
    projection(pol.choose) <- longlatproj
    pol.choose <- sp::spTransform(pol.choose,
      raster::projection(info, asText = FALSE))
  }

  #4. Calculate the number of points from each true group encompassed
  # in each calculated management unit
  if (type == "spodt") {
    info@data$estorignum <- cluster.choose@partition
    correctlabel <- rep(NA, NROW(cluster.choose@adj))
  }
  if (type == "satscan") {
    info@data$estorignum <- sp::over(info, pol.choose)$LOC_ID
    correctlabel <- rep(NA, length(pol.choose))
  }

  match.table <- tapply(info@data$estorignum, info@data$true,
    table, exclude = NULL)
  if (is.list(match.table)) {
    # Find the number of matches to each estimated partition
    match.table <- Reduce(function(...) merge(..., by = "Var1", all = TRUE),
      lapply(match.table, function(x) data.frame(x)))
    # Store the NA row for later
    remove <- which(is.na(match.table[, 1]))
    matches <- match.table
    navalues <- match.table[remove, ]
    match.table <- match.table[-remove, ]

    # Create a blank vector, with a single entry per each cluster.choose
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

    info@data$est <- factor(info@data$estorignum, labels = correctlabel)

    # Correctly label the matches table
    colnames(matches) <- c("manunit", 1:(NCOL(matches) - 1))
    rownames(matches)[1:length(correctlabel)] <- correctlabel
  } else {
    info@data$est <- 1
    correctlabel <- c("1" = 1)
    # todo: determine if this should be changed if NULL ever comes up
    matches <- NULL
  } # End if (is.list(match.table))

  #5. Calculate the amount of area each estimated polygon covers of
  # each true polygon.
  if (type == "satscan"){
    areas <- matrix(NA, nrow = length(pol.true), ncol = length(pol.choose),
      dimnames = list(paste0("true.", seq_along(pol.true)),
      paste0("choose.", seq_along(pol.choose))))
    for (ix in seq_along(pol.true)) {
      for (iy in seq_along(pol.choose)) {
        temp <- subset(pol.choose, CLUSTER == iy)
        temp <- rgeos::gIntersection(temp, pol.true[ix])
        if (is.null(temp)) {
          areas[ix, iy] <- 0
        } else {areas[ix, iy] <- rgeos::gArea(temp) / rgeos::gArea(pol.true[ix])}
        rm(temp)
      }
    }
  } else {areas <- NULL}

  #6. Plot the results if a file is specified.
    png(file.path(dir, paste0(file, ".png")),
      units = "in", res = 300, width = 8, height = 5)
    plot(info, cex = abs(min(info$omega)) + info$omega,
      pch = ifelse(info$omega >= 0, 2, 6),
      ylab = "northing",
      col = rgb(0, 0, 0, 0.4))
    if (type == "spodt"){
      for (ii in 1:length(pol.choose)) {
        lines(pol.choose[ii], lty = ii + 1, lwd = 2, xpd = TRUE)
      }
    }
    if (type == "satscan") {
      plot(pol.choose, add = TRUE, lwd = 2, xpd = TRUE)
    }
    lines(readin$lines_grouptrue, lty = 1, lwd = 0.5)
    axis(1, pos = 0)
    text(x = mean(par("usr")[1:2]), y = -400, label = "easting")
    axis(2, at = round(seq(0, extent(info)@ymax, length.out = 4) / 5, -2) * 5)
    dev.off()

  #7. Export the results.
  return(list(
    "areas" = areas,
    "matches" = matches,
    "pol.choose" = pol.choose,
    "info" = info
  ))

  #8. End the function
}
