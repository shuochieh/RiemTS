################################################################################
#### Sphere manifold: core math and S3 methods ################################
################################################################################

#' Check if x is on the sphere (within some tolerance)
#' @export
is_on_sphere = function(x, tol = 1e-6) {
  return(abs(norm(x, "2") - 1) <= tol)
}

#' Geodesic distance between two points on the sphere
#' @export
geod_sphere_core = function(x, y) {
  if (!is_on_sphere(x)) {
    stop("geod_sphere: x not on sphere")
  }
  x = x / norm(x, "2") # for numerical stability

  if (!is_on_sphere(y)) {
    stop("geod_sphere: y not on sphere")
  }
  y = y / norm(y, "2") # for numerical stability

  temp = round(c(x %*% y), 10) # resolve numerical issue (not ideal; shall revisit in the future)
  return(acos(temp))
}

#' @export
geod_core.manifold_sphere = function(mfd, x, y, ...) geod_sphere_core(x, y)

#' Exponential map on the sphere
#'
#' @param x a tangent vector
#' @param mu base point
#' @export
Exp_sphere_core = function(x, mu) {
  if (sum(abs(x)) == 0) {
    return(mu)
  }

  x_norm = sqrt(sum(x^2))
  res = cos(x_norm) * mu + sin(x_norm) * x / x_norm
  res = res / sqrt(sum(res^2)) # normalize again to avoid numerical instability

  return(res)
}

#' @export
Exp_core.manifold_sphere = function(mfd, v, mu, ...) Exp_sphere_core(v, mu)

#' Logarithmic map for the sphere
#'
#' @param x a point on the sphere
#' @param mu a base point
#' @export
Log_sphere_core = function(x, mu) {
  inner = sum(x * mu)
  inner = max(min(inner, 1), -1)
  theta = acos(inner)

  if (theta < 1e-7) {
    return(rep(0, length(mu)))
  }

  Proj = x - inner * mu
  Proj_norm = sqrt(sum(Proj^2))

  return(theta * Proj / Proj_norm)
}

#' @export
Log_core.manifold_sphere = function(mfd, x, mu, ...) Log_sphere_core(x, mu)

#' Parallel transport via geodesics for the sphere
#'
#' @param x starting point
#' @param y end point
#' @param V tangent vector at x (a vector or k by d array)
pt_sphere_core = function(x, y, V) {
  was_vector = !is.matrix(V)
  if (was_vector) {
    V = matrix(V, nrow = 1)
  }

  # 1. Check tangency for all columns simultaneously
  if (max(abs(V %*% x)) > 1e-6) {
    stop("pt_sphere: one or more column vectors in V are not tangent at x")
  }

  xy_inner = sum(x * y)

  if (xy_inner >= 1 - 1e-12) {
    if (was_vector) return(c(V))
    return(V)
  }

  e1 = x
  e2 = y - xy_inner * x
  e2 = e2 / sqrt(sum(e2^2))

  # 2. Vector projections
  a = V %*% e2

  # V_perp is the component of all k vectors orthogonal to the x-y plane.
  V_perp = V - (a %*% t(e2))
  theta = acos(xy_inner)

  # 3. Transport the basis vector e2
  e2_transported = cos(theta) * e2 - sin(theta) * e1

  # 4. Recombine
  res = (a %*% t(e2_transported)) + V_perp

  if (was_vector) {
    return(c(res))
  }

  return(res)
}

#' @export
parallel_transport.manifold_sphere = function(mfd, from, to, v, ...) {
  pt_sphere_core(from, to, v)
}

#' A basis for the tangent space at mu (sphere)
#' @export
basis_sphere = function(mu) {
  d = length(mu)
  B = matrix(svd(mu, d, d)$u[, -1], ncol = d - 1)

  return(B)
}
