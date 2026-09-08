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
geod_core.manifold_sphere = function(mfd, x, y, ...) geod_sphere(x, y, ...)

#' exponential map on the sphere
#' 
#' @param x an array of tangent vectors
#' @param mu base point (if multiple base points are supplied the number of base points should 
#'           match the number of tangent vectors)
#' 
#' @export
Exp_sphere = function (x, mu) {
  
  # Standardize inputs: x to matrix (n x q), mu to vector (q) or matrix (n x q)
  is_x_vec = is.vector(x)
  if (is_x_vec) {
    x = matrix(x, nrow = 1)
  }
  
  n = nrow(x)
  d = ncol(x)
  
  
  
  
  if (is.matrix(x)) {
    
    if (sum(abs(x)) == 0) {
      return (matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T))
    }
    
    x_norm = sqrt(rowSums(x^2))
    x_norm[which(x_norm == 0)] = 1
    std_x = x / x_norm
    res = outer(cos(x_norm), mu) + sin(x_norm) * std_x
    
    res = res / sqrt(rowSums(res^2)) # normalize again to avoid numerical instability
  } else {
    
    if (sum(abs(x)) == 0) {
      return (mu)
    }
    
    x_norm = sqrt(sum(x^2))
    res = cos(x_norm) * mu + sin(x_norm) * x / x_norm
    
    res = res / sqrt(sum(res^2))
  }
  
  return (res)
}
