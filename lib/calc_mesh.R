calc_mesh <- function(locations, boundary = prdomain,
  type = c("basic", "advanced")) {

  type <- match.arg(type, choices = c("basic", "advanced"),
    several.ok = FALSE)

  if (is.null(boundary)) {
    boundary <- INLA::inla.nonconvex.hull(locations,
      convex = -0.05)
  }

  # Create the mesh based on the type specified, where the
  # default is to use the basic mesh.
  if (type == "basic") {
    mesh <- inla.mesh.create(locations, boundary = boundary)
  }
  if (type == "advanced") {
    mesh <- inla.mesh.2d(loc = locations,
      offset = c(-0.06, -0.03),
      cutoff = 100,
      max.edge = c(100, 2000),
      boundary = boundary)
  }
  spde <- inla.spde2.matern(mesh)

  return(list("mesh" = mesh, "spde" = spde))
}
