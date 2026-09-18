library(expm)

#' Geodesic distance between two (arrays of) points on the Log-Euclidean space
#' 
#' @param x an $(m \times m)$ or $(n by m by m)$ array of SPD matrix
#' @param y an $(m \times m)$ or $(n by m by m)$ array of SPD matrix
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' y = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' geod_logE(x, y)
#' 
#' @export
geod_logE = function (x, y) {
  if (length(dim(x)) == 3) {
    nx = dim(x)[1]
  } else if (length(dim(x)) == 2) {
    nx = 1
  } else {
    stop("geod_logE: input array dimension must be 2 or 3.")
  }
  mx = dim(x)[2]

  if (length(dim(y)) == 3) {
    ny = dim(y)[1]
  } else if (length(dim(y)) == 2) {
    ny = 1
  } else {
    stop("geod_logE: input array dimension must be 2 or 3.")
  }
  my = dim(y)[2]
  
  if (nx != 1 && ny != 1 && nx != ny) {
    stop("geod_logE: number of input matrices must match or equal to 1; otherwise broadcasting is undefined.")
  } 
  if (mx != my) {
    stop("geod_logE: input matrix dimensions must match.")
  }
  
  m = mx
  
  if (nx == 1) {
    log_x = logm(x)
  } else {
    log_x = array(NA, dim = c(nx, m, m))
    for (i in 1:nx) {
      log_x[i,,] = logm(x[i,,])
    }
  }
  
  if (ny == 1) {
    log_y = logm(y)
  } else {
    log_y = array(NA, dim = c(ny, m, m))
    for (i in 1:ny) {
      log_y[i,,] = logm(y[i,,])
    }
  }
  
  if (nx == 1 && ny == 1) {
    D = log_x - log_y
    return (sqrt(sum(D^2)))
  } else if (nx == 1) {
    vec_x = c(log_x)
    mat_y = matrix(log_y, nrow = ny, ncol = m * m)
    D = t(t(mat_y) - vec_x)
    return (sqrt(rowSums(D^2)))
  } else if (ny == 1) {
    mat_x = matrix(log_x, nrow = nx, ncol = n * n)
    vec_y = c(log_y)
    D = t(t(mat_x) - vec_y)
    return (sqrt(rowSums(D^2)))
  } else {
    mat_x = matrix(log_x, nrow = nx, ncol = m * m)
    mat_y = matrix(log_y, nrow = ny, ncol = m * m)
    return (sqrt(rowSums((mat_x - mat_y)^2)))
  }
}

#' @export
geod.manifold_logEuclidean = function (mfd, x, y) {
  geod_logE(x, y)
}

#' helper function to compute the first difference ratio
#' 
#' @param a a vector of real numbers (length > 1)
#' @param type either "exp" or "log"
#' @param tol threshold for numerical approximation
#' 
first_diff_ratio = function (a, type = "exp", tol = 1e-5) {
  L = length(a)
  if (L <= 1) {
    stop("first_diff_ratio: argument must have length > 1")
  }
  
  Delta = outer(a, a, "-")
  near_equal = abs(Delta) < tol
  
  if (type == "exp") {
    ea = exp(a)
    ea_mat = outer(ea, ea, "-")
    res = ea_mat / Delta
    
    ea_j = matrix(ea, nrow = L, ncol = L, byrow = TRUE)
    taylor = ea_j * (1 + Delta / 2 + (Delta^2) / 6 + (Delta^3) / 24)
    
    res[near_equal] = taylor[near_equal]
  } else if (type == "log") {
    if (any(a <= 0)) {
      stop("first_diff_ratio: log requires strictly positive values.")
    }
    la = log(a)
    la_mat = outer(la, la, "-")
    res = la_mat / Delta
    
    a_j = matrix(a, nrow = L, ncol = L, byrow = TRUE)
    h = Delta / a_j
    taylor = (1 / a_j) * (1 - h / 2 + (h^2) / 3 - (h^3) / 4)
    
    res[near_equal] = taylor[near_equal]
  }
  
  return (res)
}

#' Differentials of the matrix exp and log
#' Computes Dexp_P(Q) or Dlog_P(Q)
#' 
#' @param P an $m \times m$ (SPD) matrix (base)
#' @param Q an $m \times m$ matrix
#' @param type either "exp" or "log"
#' 
#' @export
diff_explog = function (P, Q, type = "exp") {
  
  spec_decomp = eigen(P, symmetric = TRUE)
  U = spec_decomp$vectors
  lambda = spec_decomp$values
  
  Q_tilde = t(U) %*% Q %*% U
  FD_mat = first_diff_ratio(lambda, type = type)
  
  res = U %*% (FD_mat * Q_tilde) %*% t(U)
  
  res = 0.5 * (res + t(res))
  
  return (res)
}

#' exponential map on the log-Euclidean geometry
#' 
#' @param z an $m \times m$  or $n \times m \times m$ array of tangent vectors (identified as symmetric matrices)
#' @param x an $m \times m$ SPD matrix (as the base point)
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' z = crossprod(matrix(rnorm(9), 3, 3)) * 0.1
#' Exp_logE(z, x)
#' 
#' @export
Exp_logE = function (z, x) {
  if (length(dim(z)) == 3) {
    n = dim(z)[1]
    was_matrix = FALSE
  } else if (length(dim(z)) == 2) {
    n = 1
    z = array(z, dim = c(1, dim(z)))
    was_matrix = TRUE
  } else {
    stop("Exp_logE: number of dimensions of z must be either 2 or 3.")
  }
  
  m = dim(z)[2]
  if (dim(z)[3] !=m || !all(dim(x) == c(m, m))) {
    stop("Exp_logE: input matrix dimension must match (m x m).")
  }
  
  res = array(NA, dim = c(n, m, m))
  log_x = logm(x)
  for (i in 1:n) {
    res[i,,] = expm(log_x + diff_explog(x, z[i,,], "log"))
  }
  
  if (was_matrix) {
    res = res[1,,]
  }
  return (res)
}

#' @export
Exp_mfd.manifold_logEuclidean = function (mfd, p, v, ...) {
  Exp_logE(z = v, x = p)
}

#' logarithm map on the log-Euclidean geometry
#' 
#' @param x an $m \times m$ SPD matrix (as the base point)
#' @param y an $m \times m$  or $n \times m \times m$ array of SPD matrices
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' y = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' Log_logE(x, y)
#'  
#' @export
Log_logE = function (x, y) {
  if (length(dim(y)) == 3) {
    n = dim(y)[1]
    was_matrix = FALSE
  } else if (length(dim(y)) == 2) {
    n = 1
    y = array(y, dim = c(1, dim(y)))
    was_matrix = TRUE
  } else {
    stop("Log_logE: number of dimensions of y must be either 2 or 3.")
  }
  
  m = dim(y)[2]
  if (dim(y)[3] != m || !all(dim(x) == c(m, m))) {
    stop("Log_logE: input matrix dimension must match (m x m).")
  }
  
  res = array(NA, dim = c(n, m, m))
  log_x = logm(x)
  for (i in 1:n) {
    res[i,,] = diff_explog(log_x, logm(y[i,,]) - log_x, "exp")
  }
  if (was_matrix) {
    res = res[1,,]
  }
  
  return (res)
}

#' @export
Log_mfd.manifold_logEuclidean = function (mfd, p, q, ...) {
  Log_logE(x = p, y = q)
}






