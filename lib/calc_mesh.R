calc_mesh <- function(locations, boundary = prdomain) {
  mesh <- inla.mesh.2d(loc = locations,
    offset = c(-0.06, -0.03),
    cutoff = 100,
    max.edge = c(100, 2000),
    boundary = boundary)

  spde <- inla.spde2.matern(mesh)

  return(list("mesh" = mesh, "spde" = spde))
}
