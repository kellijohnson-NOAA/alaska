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
#' @param projection The spatial projection of the results. The default is
#' \code{sp::CRS("+proj=longlat +ellps=WGS84")}, but any projection can
#' be used. If latitude and longitude is used the function will create a
#' shapefile for the results. If cartesian coordinates are used, then the
#' results will be read in and a hull is created around the points that
#' are selected for each group because SaTScan does not create a shapefile
#' for cartesian coordinates.
#' @param criticalvalue A numeric value indicating the threshold limit for the
#' pvalue, where above the threshold the estimated polygons will be eliminated
#' from the results. The default value is 0.05.
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
  projection = sp::CRS("+proj=longlat +ellps=WGS84"),
  criticalvalue = 0.05, ...) {
  type <- "satscan"
  if (is.null(exe)) type <- "spodt"

  #1. Pull out the relevant portions of readin
  # info is a subset of the model results for Omega,
  # that specify the point locations and the estimated values
  # that lie within the inner boundary of the mesh.
  # These were determined using \code\link{{calc_meshbound}}
  info <- readin$info

  #2. Calculate the true polygons
  bounds <- calc_meshbound(readin$mesh,
    projection = sp::CRS(raster::projection(readin$info)))
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
  }
  if (type == "satscan") {
    export_prm(file_in = prmtemplate, dir = dir, file = file, ...)
    export_geo(points = info, dir = dir, file = file,
      projection = projection)
    # Double check the direction of slashes and call system
    system(paste(
      gsub("/", "\\\\", exe),
      file.path(gsub("/", "\\\\", dir), paste0(file, ".prm"), fsep = "\\")
    ), show.output.on.console = FALSE)
    # Read in the results
    textfile <- read_txt(dir = dir, file = file)
    # Get p values
    if (is.data.frame(textfile$coordinates) | is.matrix(textfile$coordinates)) {
      pvalues <- as.numeric(sapply(strsplit(
        textfile$coordinates[grep("P-", textfile$coordinates[, 1]), ],
        ":"), "[[", 2))
    }
    # If the shapefile exists, read it in.
    shapefile <- file.path(dir, paste0(file, ".col"))
  }

  # split the full polygon by the estimated polygons
  if (exists("lines.choose")) {
    if (length(unique(cluster.choose@partition)) == 1) {
      pol.choose <- bounds$hull
    } else {
      pol.choose <- calc_polys(rgeos::gBuffer(bounds$hull, width = 1),
        lines.choose)
    }
  } else {
      pol.choose <- SpatialPolygons(
        lapply(textfile$cluster,
        function(x, data = info@coords) {
          X <- data[x, ]
          X.chull <- chull(X)
          X.chull <- c(X.chull, X.chull[1])
          Polygons(list(Polygon(X[X.chull, ])), parent.frame()$i[])
      }))
    projection(pol.choose) <- raster::projection(info, asText = FALSE)
    if (file.exists(paste0(shapefile, ".shp"))) {
      pol.hull <- readShapePoly(shapefile)
      pol.hull <- SpatialPolygons(pol.hull@polygons)
      projection(pol.hull) <- raster::projection("+proj=longlat +ellps=WGS84")
      if (!exists("pol.choose")) pol.choose <- pol.hull
    }
    if (exists("pvalues")) {
      pol.choose <- pol.choose[pvalues <= criticalvalue]
    }
  }

  #4. Calculate the number of points from each true group encompassed
  # in each calculated management unit
  if (length(pol.choose) > 0) {
    temp <- sp::over(info, pol.choose)
    correctlabel <- rep(NA, length(pol.choose))
    if ("LOC_ID" %in% colnames(temp)) {
        info@data$estorignum <- temp$LOC_ID
    } else {info@data$estorignum <- temp}
    rm(temp)
  } else {
    info@data$estorignum <- NA
    info@data$est <- NA
  }

  match.table <- tapply(info@data$estorignum, info@data$true,
    table, exclude = NULL)
  if (is.list(match.table)) {
    # Find the number of matches to each estimated partition
    match.table <- Reduce(function(...) merge(..., by = "Var1", all = TRUE),
      lapply(match.table, function(x) data.frame(x)))
    # Store the NA row for later
    matches <- match.table
    match.table <- match.table[-which(is.na(match.table[, "Var1"])), ]

    # Create a blank vector, with a single entry per each cluster.choose
    names(correctlabel) <- match.table[, 1]
    rownames(match.table) <- match.table[, 1]

   # Remove the names and rename the columns
    match.table[, "Var1"] <- 0
    colnames(match.table) <- 0:(NCOL(match.table) - 1)

    while (NCOL(match.table) > 1) {
      if (max(match.table, na.rm = TRUE) == 0) {break}
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
    if (NCOL(matches) > 2) {
      matches <- reshape(matches, direction = "long",
        idvar = "Var1",
        varying = grep("Freq", colnames(matches), value = TRUE))
    } else {
      matches <- data.frame("Var1" = matches$Var1, "time" = "x",
        "Freq" = matches$Freq)
    }
    if (!all(matches$time %in% c("x", "y"))) {
      stop("There is an additional variable name in", matches$time)
    } else {
      matches$time <- as.numeric(factor(matches$time))
    }
    matches <- merge(x = matches, y = data.frame(table(info@data$true)),
      by.x = "time", by.y = "Var1", all = TRUE)
    colnames(matches)[which(colnames(matches) == "Freq.y")] <- "true"
    colnames(matches) <- gsub("\\.x", "", colnames(matches))
    matches$Var1 <- factor(matches$Var1, labels = correctlabel)
  } else {
    matches <- data.frame(
      "time" = as.numeric(as.character(names(match.table))),
      "Var1" = NA, "Freq" = match.table, "true" = sum(match.table))
  } # End if (is.list(match.table))
  matches$replicate <- readin$replicate
  matches$scenario <- gsub("sim_|\\.RData", "", basename(readin$file))
  matches$subpopulations <- strsplit(matches$scenario, "-")[[1]][1]
  matches$alpha <- strsplit(matches$scenario, "-")[[1]][3]
  matches$Freq[is.na(matches$Freq)] <- 0
  matches$decimal <- matches$Freq / matches$true
  matches$decimal <- ifelse(matches$Var1 == matches$time |
    is.na(matches$Var1), 1, -1) * matches$decimal
  matches$decimal <- ifelse(is.na(matches$Var1), -1, 1) * matches$decimal

  #5. Calculate the amount of area each estimated polygon covers of
  # each true polygon.
  if (type == "satscan" & length(pol.choose) > 0){
    # Use the hull polygon rather than the outer bounds
    truth <- bounds$hull
    if (length(pol.true) == 1) {
       truth <- bounds$hull
    } else {truth <- calc_polys(bounds$hull, readin$lines_grouptrue)}

    areas <- matrix(NA, nrow = length(truth), ncol = length(pol.choose),
      dimnames = list(paste0("true.", seq_along(truth)),
      paste0("choose.", seq_along(pol.choose))))
    for (ix in seq_along(truth)) {
      for (iy in seq_along(pol.choose)) {
        if (class(pol.choose) == "SpatialPolygons") {
          temp <- pol.choose[iy]
        } else {temp <- subset(pol.choose, CLUSTER == iy)}
        if (.hasSlot(temp@polygons[[1]]@Polygons[[1]], "area")) {
          if (temp@polygons[[1]]@Polygons[[1]]@area == 0) {
          temp <- NULL
          }} else {
            temp <- rgeos::gIntersection(temp, truth[ix])
          } #End if area == 0
        if (is.null(temp)) {
          areas[ix, iy] <- 0
        } else {areas[ix, iy] <- rgeos::gArea(temp) / rgeos::gArea(truth[ix])}
        rm(temp)
      }
    }
  } else {areas <- NULL}

  #6. Plot the results if a file is specified.
    png(file.path(dir, paste0(file, ".png")),
      units = "in", res = 300, width = 8, height = 5)
    par(xpd = TRUE)
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
    text(x = mean(par("usr")[1:2]), y = -350, label = "easting")
    axis(2, at = round(seq(0, extent(info)@ymax, length.out = 4) / 5, -2) * 5)
    dev.off()

  #7. Export the results.
  return(list(
    "areas" = areas,
    "matches" = matches,
    "pol.choose" = pol.choose,
    "info" = info,
    "n_years" = readin$n_years
  ))

  #8. End the function
}
