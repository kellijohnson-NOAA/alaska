###############################################################################
###############################################################################
## Purpose:    Stock analysis of cod in Alaska
##             Create mesh from INLA
## Author:     Kelli Faye Johnson
## Contact:    kellifayejohnson@gmail.com
## Date:       2016-08-23
## Comments:
###############################################################################
###############################################################################

###############################################################################
# Create the mesh using the projected data so all measurements are in m
# The function inla.mesh.2d() is recommended and will create
# Constrained Refined Delaunay Triangulation (CRDT) == mesh
# The boundary will have a variance 2x the variance of the domain
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
###############################################################################

# use a mesh based on a non-convex hull
# to avoid adding many small triangles outside the domain of interest
# (more triangles = larger computation times)
# inla.nonconvex.hull(points, convex, concave, resolution, eps)
# points = 2D point coordinates (2-column matrix)
# convex = desired extension radius
#   Determines the smallest allowed convex curvature radius.
#   Negative values are fractions of the approximate initial set diameter.
# concave = minimal concave curvature radius
#   Default is concave = convex.
# resolution = internal computation resolution
#   A warning will be issued when this needs to be increased for
#   higher accuracy, with the required resolution stated.
# eps = The polygonal curve simplification tolerance used for simplifying
#   the resulting boundary curve.
prdomain <- inla.nonconvex.hull(coordinates(data.all),
  convex = -0.05, resolution = c(40, 15))

# inla.mesh.2d(loc = NULL, loc.domain = NULL, offset = NULL,
#   n = NULL, boundary = NULL, interior = NULL, max.edge,
#   min.angle = NULL, cutoff = 1e-12, plot.delay = NULL)
# loc = locations used as initial triangulation nodes
# loc.domain = a polygon to determine the domain extent
# max.edge = specifies the maximum allowed triangle edge lengths
#   in the inner domain and in the outer extension.
#   A scalar or length two vector, on the SAME SCALE UNIT as
#   the coordinates.
# offset = a numeric, or length two vector.
#   If negative then it is a factor relative to the approx. data diameter.
#   If positive it is the extension distance on same scale as coords.
# n = initial number of points on the extended boundary
# interior = a list of segments to specify interior constraints,
#   each one of inla.mesh.segment class.
#   A good mesh needs to have triangles as regular as possible
#   in size and shape.
# min.angle = a scalar or length two vector,
#   to specify the minimum internal angles of the triangles on
#   the inner domain and on the outer extension.
#   Values up to 21 guarantee the convergence of the algorithm.
# cutoff = minimum allowed distance between points.
#   Points at a closer distance than the supplied value are replaced
#   by a single vertex to avoid small triangles.
#   Critical when we have some very close points, either for
#   point locations or on the domain boundary.
mesh <- inla.mesh.2d(loc = coordinates(data.all),
  offset = c(-0.06, -0.03),
  cutoff = 100,
  max.edge = c(100, 2000),
  boundary = prdomain)

spde <- inla.spde2.matern(mesh)

# End of file
