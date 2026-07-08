################################################################################
#### Generic parallel transport ################################################
################################################################################
#
# Gives pt_bws/pt_gr_core/pt_sphere a uniform calling convention. The
# underlying math is genuinely different per manifold (BWS integrates an
# ODE along the geodesic; Grassmannian and sphere use closed-form
# single-step formulas), so only the *interface* is unified here - each
# method's implementation lives alongside that manifold's other core math
# (R/bws.R, R/grassmann.R, R/sphere.R).

#' Parallel transport of a tangent vector along a geodesic
#'
#' @param mfd a manifold object
#' @param from starting point
#' @param to end point
#' @param v tangent vector at `from`
#' @param ... manifold-specific extra arguments (e.g. `method` for manifold_bws())
#' @return the transported tangent vector, at `to`
#' @examples
#' mfd = manifold_sphere()
#' x = c(1, 0, 0)
#' y = c(0, 1, 0)
#' v = c(0, 0, 1)   # tangent to the sphere at x
#' parallel_transport(mfd, x, y, v)
#' @export
parallel_transport = function(mfd, from, to, v, ...) UseMethod("parallel_transport")
