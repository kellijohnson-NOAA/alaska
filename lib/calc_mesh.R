calc_mesh <- function(locations, boundary = prdomain,
  type = c("basic", "advanced", "default")) {

  type <- match.arg(type, choices = c("basic", "advanced", "default"),
    several.ok = FALSE)

  if (is.null(boundary)) {
    boundary <- INLA::inla.nonconvex.hull(locations,
      convex = -0.05)
  }

  # Find the area of the bounding box and take
  # 5% of the square root of the area
  # as the cutoff for the created mesh
  if (is.data.frame(locations)) {
    cutoff <- t(sp::bbox(array(locations)))
  } else {
    cutoff <- t(sp::bbox(locations))
  }
  cutoff <- rbind(
    c(cutoff[1, 1], cutoff[2, 2]),
    cutoff[1, ],
    c(cutoff[2, 1], cutoff[1, 2]),
    cutoff[2, ])
  rownames(cutoff) <- NULL
  cutoff <- SpatialPolygons(list(Polygons(list(Polygon(list(cutoff))), ID = 1)))
  cutoff <- cutoff@polygons[[1]]@area
  cutoff <- sqrt(ceiling(cutoff)) * 0.05

  # Create the mesh based on the type specified, where the
  # default is to use the basic mesh.
  if (type == "basic") {
    mesh <- inla.mesh.create(locations, boundary = boundary,
      cutoff = cutoff)
  }
  if (type == "advanced") {
    mesh <- inla.mesh.2d(loc = locations,
      offset = c(-0.06, -0.03),
      cutoff = 100,
      max.edge = c(100, 2000),
      boundary = boundary)
  }
  if (type == "default") {
    mesh <- inla.mesh.create(locations)
  }

  return(list("mesh" = mesh, "cutoff" = cutoff))
}
