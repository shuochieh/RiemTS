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
Exp_mfd.manifold_sphere = function(mfd, v, mu, ...) Exp_sphere(v, mu, ...)

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
Log_mfd.manifold_sphere = function(mfd, x, mu, ...) Log_sphere(x, mu, ...)





