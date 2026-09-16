library(Matrix)

#' Check if x is on the sphere (within some tolerance)
#' @param x a numeric vector or matrix
#' @param tol tolerance on `|| x || - 1`
#' @return TRUE/FALSE
#' @examples
#' is_on_sphere(c(1, 0, 0))
#' is_on_sphere(c(1, 1, 0))
#' @export
is_on_sphere = function(x, tol = 1e-6) {
  # Treat vector as 1-row matrix to unify norm calculation
  if (is.vector(x)) {
    x = matrix(x, nrow = 1)
  }
  
  if (!is.matrix(x)) {
    stop("is_on_sphere: input must be a vector or matrix.")
  }
  
  # Compute Euclidean norms along rows
  norms = sqrt(rowSums(x^2))
  
  # Returns logical scalar (for vector) or logical vector (for matrix)
  res = abs(norms - 1) <= tol
  if (length(res) == 1) {
    return(res[1])
  }
  
  return(res)
}

#' Geodesic distance between two (array) of points on the sphere
#' If x and y are both arrays, a vector of n distances corresponding to pointwise
#' geodesic distances is returned.
#' 
#' @param x an $(n \times q)$ array of points or a $q$-dimensional vector
#' @param y an $(n \times q)$ array of points or a $q$-dimensional vector
#' @tol tolerance in checking sphere membership
#' 
#' @examples 
#' x = matrix(c(1,0,0,1,0,0), ncol = 3)
#' y = c(0, 0, 1)
#' geod_sphere(x, y)
#' 
#' @export
geod_sphere = function (x, y, tol = 1e-6) {
  
  # Check sphere membership 
  if (!all(is_on_sphere(x, tol = tol))) {
    stop("geod_sphere: x contains points not on sphere")
  }
  if (!all(is_on_sphere(y, tol = tol))) {
    stop("geod_sphere: y contains points not on sphere")
  }
  
  # standardize input to matrices
  if (is.vector(x)) {
    x = matrix(x, nrow = 1)
  }
  if (is.vector(y)) {
    y = matrix(y, nrow = 1)
  }
  
  if (!is.matrix(x) || !is.matrix(y)) {
    stop("geod_sphere: x and y must be vectors or 2D matrices.")
  }
  
  # handle broadcasting
  nx = nrow(x)
  ny = nrow(y)
  if (nx == 1 && ny > 1) {
    x = matrix(x, nrow = ny, ncol = ncol(x), byrow = TRUE)
  }
  if (ny == 1 && nx > 1) {
    y = matrix(y, nrow = nx, ncol = ncol(y), byrow = TRUE)
  }
  
  if (nrow(x) != nrow(y) || ncol(x) != ncol(y)) {
    stop("geod_sphere: dimension mismatch between x and y.")
  }
  
  # Compute row norms efficiently
  norm_x = sqrt(rowSums(x^2))
  norm_y = sqrt(rowSums(y^2))
  
  # Normalize to enforce exact unit length for safety
  x = x / norm_x
  y = y / norm_y
  
  # compute the geodesic distance
  dot_prod = rowSums(x * y)
  dot_prod = pmax(pmin(dot_prod, 1), -1) # Clamp values to [-1, 1] to avoid NaN from floating-point overshoot
  res = acos(dot_prod)
  
  return (res)
}

#' @export
geod.manifold_sphere = function(mfd, x, y, ...) geod_sphere(x, y, ...)

#' exponential map on the sphere
#' 
#' @param x an array of tangent vectors
#' @param mu base point (if multiple base points are supplied the number of base points should 
#'           match the number of tangent vectors)
#' @param tol tolerance parameter for numerical approximation near zero
#' 
#' @examples 
#' mu = c(0, 1)
#' x = c(pi/2, 0)
#' Exp_sphere(x, mu)
#' 
#' @export
Exp_sphere = function (x, mu, tol = 1e-4) {
  
  # Standardize inputs: x to matrix (n x q), mu to vector (q) or matrix (n x q)
  is_x_vec = is.vector(x)
  if (is_x_vec) {
    x = matrix(x, nrow = 1)
  }
  
  n = nrow(x)
  d = ncol(x)
  
  if (is.vector(mu)) {
    if(length(mu) != d) {
      stop("Exp_sphere: dimension mismatch between x and mu.")
    }
    
    mu_mat = matrix(mu, nrow = n, ncol = d, byrow = TRUE)
  } else if (is.matrix(mu)) {
    if (nrow(mu) != n || ncol(mu) != d) {
      stop("Exp_sphere: dimension mismatch between x and mu.")
    }
    
    mu_mat = mu
  } else {
    stop("Exp_sphere: mu must be a vector or a matrix.")
  }
  
  # Squared norm avoids unnecessary sqrt for small vectors
  r2 = rowSums(x^2)
  r  = sqrt(r2)
  
  # Taylor expansion threshold at tol
  # Below tol, 1 - r2/6 + r2^2/120 approximates sin(r)/r 
  small_mask = (r < tol)
  
  sinc = numeric(n)
  
  # Standard division for normal/large vectors
  sinc[!small_mask] = sin(r[!small_mask]) / r[!small_mask]
  
  # Smooth, highly accurate Taylor series for small/zero vectors
  sinc[small_mask]  = 1 - r2[small_mask] / 6 + (r2[small_mask]^2) / 120
  
  # Pointwise Exponential Map
  res = cos(r) * mu_mat + sinc * x
  
  # Safeguard unit-norm mapping
  res = res / sqrt(rowSums(res^2))
  
  if (is_x_vec) {
    return(res[1, ])
  }
  
  return(res)
}

#' @export
Exp_mfd.manifold_sphere = function (mfd, v, mu, ...) Exp_sphere(v, mu, ...)

#' logarithmic map for the sphere
#' 
#' @param x an array of points on the sphere
#' @param mu base point(s)
#' @param tol tolerance parameter for numerical approximation near zero
#' @param tol_antipodal tolerance parameter for checking antipodal points 
#' 
#' @examples 
#' p1 = c(1, 0, 0)
#' p2 = c(0, 1, 0)
#' p3 = c(0, 0, 1)
#' x = rbind(p1, p2, p3)
#' mu = rbind(p3, p2, p1)
#' Log_sphere(x, mu)
#' 
#' @export
Log_sphere = function (x, mu, tol = 1e-4, tol_antipodal = 1e-7) {
  
  x_was_vec = is.vector(x)
  mu_was_vec = is.vector(mu)
  
  if (x_was_vec) {
    x = matrix(x, nrow = 1)
  }
  if (mu_was_vec) {
    mu = matrix(mu, nrow = 1)
  }
  
  if (!is.matrix(x) || !is.matrix(mu)) {
    stop("Log_sphere: x and mu must be vectors or matrices.")
  }
  
  nx = nrow(x)
  nmu = nrow(mu)
  d = ncol(x)
  
  if (ncol(mu) != d) {
    stop("Log_sphere: dimension mismatch (columns of x and mu).")
  }
  
  if (nmu == 1 && nx > 1) {
    mu = matrix(mu, nrow = nx, ncol = d, byrow = TRUE)
    n = nx
  } else if (nx == 1 && nmu > 1) {
    x = matrix(x, nrow = nmu, ncol = d, byrow = TRUE)
    n = nmu
  } else if (nx == nmu) {
    n = nx
  } else {
    stop("Log_sphere: row dimensions must match or one input must have 1 row.")
  }
  
  # inner product
  inner = rowSums(x * mu)
  inner = pmax(pmin(inner, 1), -1)
  
  theta = acos(inner)
  
  if (any(abs(theta - pi) < tol)) {
    warning("Log_sphere:contains antipodal pairs where Log map is undefined.")
  }
  
  Proj = x - inner * mu
  
  scale = numeric(n)
  small_mask = (theta < tol)
  
  scale[!small_mask] = theta[!small_mask] / sin(theta[!small_mask])
  
  theta2 = theta[small_mask]^2
  scale[small_mask] = 1 + theta2 / 6 + (7 * theta2^2) / 360
  
  res = scale * Proj
  
  if (x_was_vec && mu_was_vec) {
    return (c(res))
  }
  
  return (res)
}

#' @export
Log_mfd.manifold_sphere = function (mfd, x, mu, ...) Log_sphere(x, mu, ...)

#' Tangency check for the sphere
#' 
tangency_check_sphere = function (x, V, tol = 1e-8) {
  if (is.vector(V)) {
    V = matrix(V, nrow = 1)
  }
  
  if (!is.null(tol) && !is.infinite(tol)) {
    inner_Vx = c(V %*% x)
    if (any(sqrt(sum(inner_Vx^2)) > tol)) {
      return (FALSE)
    }
  }
  return (TRUE)
}

#' Parallel transport via geodesics for the sphere
#'
#' @param x starting point
#' @param y end point
#' @param V tangent vector at x (a vector or k by d array)
#' @param tol tolerance for tangency check (default is 1e-8)
#' 
#' @examples 
#' pt_sphere(c(0,0,1), c(0,sin(pi/5),cos(pi/5)), c(1,-1,0))
#' 
#' @export
pt_sphere = function(x, y, V, tol = 1e-8) {
  was_vector = is.vector(V)
  if (was_vector) {
    V = matrix(V, nrow = 1)
  }
  
  k = nrow(V)
  d = ncol(V)
  
  if (length(x) != d || length(y) != d) {
    stop("pt_sphere: x, y, and columns of V must have matching dimension.")
  }
  
  # Tangency check
  if (!tangency_check_sphere(x, V, tol = tol)) {
    stop("pt_sphere: one or more row vectors in V are not tangent at x")
  }
  
  xy_inner = sum(x * y)
  xy_inner = pmax(pmin(xy_inner, 1), -1)
  
  if (xy_inner >= 1 - 1e-12) {
    # identity transport
    if (was_vector) {
      return (c(V)) 
    }
    return (V)
  }
  
  if (xy_inner <= -1 + 1e-12) {
    stop("pt_sphere: x and y are (nearly) antipodal; parallel transport along unique geodesic is not defined.")
  }
  
  theta = acos(xy_inner)
  
  e1 = x
  proj_y = y - xy_inner * x
  proj_y_norm = sqrt(sum(proj_y^2))
  e2 = proj_y / proj_y_norm
  
  a = V %*% e2
  V_perp = V - a %*% matrix(e2, nrow = 1)
  
  e2_transported = cos(theta) * e2 - sin(theta) * e1
  res = (a %*% matrix(e2_transported, nrow = 1)) + V_perp
  
  if (was_vector) {
    return (c(res))
  }
  return (res)
}

#' @export
ptransport.manifold_sphere = function(mfd, from, to, v, ...) {
  pt_sphere(from, to, v, ...)
}

#' A basis for the tangent space at mu (sphere)
#' @param mu base point on the sphere
#' @return a (d by d-1) matrix whose columns form an orthonormal basis of the tangent space at `mu`
#' 
#' @examples
#' basis_sphere(c(1, 0, 0))
#' 
#' @export
basis_sphere = function(mu) {
  d = length(mu)
  B = matrix(svd(mu, d, 1)$u[, -1], ncol = d - 1)
  
  return(B)
}

#' @export
basis.manifold_sphere = function (mfd, mu) {
  basis_sphere(mu)
}

#' Compute the Riemannian Hessian vector action H[v] on the sphere
#' 
#' Evaluates the action of the Riemannian Hessian of f(x) = 0.5 * d^2(x, mu)
#' on one or more tangent vectors v
#' 
#' @param x base point where the Hessian is evaluated
#' @param mu target reference point
#' @param V tangent vector(s) at x (vector of length d or k by d matrix)
#' @param tol tolerance for tangency check
#' 
#' @export
Hess_sphere = function (x, mu, V, tol = 1e-8) {
  was_vector = is.vector(V)
  if (was_vector) {
    V = matrix(V, nrow = 1)
  }
  
  k = nrow(V)
  d = ncol(V)
  
  if (length(x) != d || length(mu) != d) {
    stop("Hess_sphere: x, mu, and columns of V must have matching dimension")
  }
  
  # Tangency check
  if (!tangency_check_sphere(x, V, tol = tol)) {
    stop("pt_sphere: one or more row vectors in V are not tangent at x")
  }
  
  cos_theta = sum(x * mu)
  cos_theta = pmax(pmin(cos_theta, 1), -1)
  theta = acos(cos_theta)
  theta2 = theta^2
  
  if (theta < 1e-4) {
    c2 = theta2 / 3 + (theta2^2) / 45
    c1 = 1 - c2
  } else {
    c1 = theta / tan(theta)
    c2 = 1 - c1
  }
  
  proj_x = x - cos_theta * mu
  sin_theta = sqrt(sum(proj_x^2))
  if (sin_theta < 1e-12) {
    if (was_vector) {
      return (c(V)) # Hessian is identity around mu
    } 
    return (V)
  }
  
  u = proj_x / sin_theta
  
  v_mu = c(V %*% mu)
  v_u = c(V %*% u)
  
  term1 = c1 * (V - v_mu %*% matrix(mu, nrow = 1))
  term2 = c2 * (v_u %*% matrix(u, nrow = 1))
  
  res = term1 + term2
  if (was_vector) {
    return(c(res))
  }
  return (res)
}

#' @export
Hessian.manifold_sphere = function (mfd, x, mu, V, ...) {
  Hess_sphere(x, mu, V, ...)
}




