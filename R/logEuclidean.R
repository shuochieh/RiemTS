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
    mat_x = matrix(log_x, nrow = nx, ncol = m * m)
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
geod.manifold_logE = function (mfd, x, y, ...) {
  geod_logE(x, y)
}

#' helper function to compute the first difference ratio
#' 
#' @param a a vector of real numbers (length > 1)
#' @param type either "exp" or "log"
#' @param tol threshold for numerical approximation
#' 
#' @noRd
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

#' Exponential map on the log-Euclidean geometry
#' 
#' @param z an $m \times m$  or $n \times m \times m$ array of tangent vectors (identified as symmetric matrices)
#' @param x an $m \times m$ SPD matrix (as the base point)
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' z = crossprod(matrix(rnorm(9), 3, 3)) * 0.1
#' Exp_logE(z, x)
#' 
#' @references 
#'   Arsigny, V., Fillard, P., Pennec, X., & Ayache, N. (2007). 
#'   Geometric means in a novel vector space structure on symmetric positive-definite matrices. 
#'   \emph{SIAM Journal on Matrix Analysis and Applications}, 29(1), 328--347. 
#'   
#'   Wagner vom Berg, G., Röhr, V., Platt, D., & Blankertz, B. (2024). 
#'   A new canonical log-Euclidean kernel for symmetric positive definite matrices for EEG analysis. 
#'   \emph{IEEE Transactions on Biomedical Engineering}. 
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

#' @rdname Exp_mfd
#' @export
Exp_mfd.manifold_logE = function (mfd, p, v, ...) {
  Exp_logE(z = v, x = p)
}

#' Logarithm map on the log-Euclidean geometry
#' 
#' @param x an \eqn{m \times m} SPD matrix (as the base point)
#' @param y an \eqn{m \times m}  or \eqn{n \times m \times m} array of SPD matrices
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' y = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' Log_logE(x, y)
#'  
#' @references 
#'   Arsigny, V., Fillard, P., Pennec, X., & Ayache, N. (2007). 
#'   Geometric means in a novel vector space structure on symmetric positive-definite matrices. 
#'   \emph{SIAM Journal on Matrix Analysis and Applications}, 29(1), 328--347. 
#'   
#'   Wagner vom Berg, G., Röhr, V., Platt, D., & Blankertz, B. (2024). 
#'   A new canonical log-Euclidean kernel for symmetric positive definite matrices for EEG analysis. 
#'   \emph{IEEE Transactions on Biomedical Engineering}. 
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
Log_mfd.manifold_logE = function (mfd, p, q, ...) {
  Log_logE(x = p, y = q)
}

#' Parallel transport along geodesic on the log-Euclidean geometry
#' 
#' @param p an \eqn{m \times m} SPD matrix (start)
#' @param q an \eqn{m \times m} SPD matrix (end)
#' @param x an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric matrices, 
#'          identified as tangent vectors at p
#' 
#' @examples 
#' p = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' q = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' v = matrix(rnorm(9), ncol = 3)
#' v = 0.5 * (v + t(v))
#' pt_logE(p, q, v)
#' 
#' @references Thanwerdas, Y., and Pennec, X. (2023). O(n)-invariant Riemannian 
#'   metrics on SPD matrices. \emph{Linear Algebra and its Applications}, 
#'   661: 163-201. 
#' 
#' @export
pt_logE = function (p, q, x) {
  if (length(dim(x)) == 3) {
    n = dim(x)[1]
    was_matrix = FALSE
  } else if (length(dim(x)) == 2) {
    n = 1
    x = array(x, dim = c(1, dim(x)))
    was_matrix = TRUE
  } else {
    stop("pt_logE: number of dimensions of x must be either 2 or 3.")
  }
  
  m = dim(x)[2]
  res = array(NA, dim = c(n, m, m))
  
  log_q = logm(q)
  for (i in 1:n) {
    temp = diff_explog(P = p, Q = x[i,,], type = "log")
    res[i,,] = diff_explog(P = log_q, Q = temp, type = "exp")
  }
  if (was_matrix) {
    return (res[1,,])
  }
  return (res)
}

#' @export
ptransport.manifold_logE = function (mfd, from, to, v, ...) {
  pt_logE(p = from, q = to, x = v)
}

#' A basis for the tangent space at p (SPD matrix in the Log-Euclidean geometry)
#'
#' @param p an \eqn{m \times m} SPD matrix (base point)
#' 
#' @examples 
#' p = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' basis_logE(p)
#' 
#' @export
basis_logE = function (p) {
  m = dim(p)[1]
  if (length(dim(p)) != 2 || dim(p)[2] != m) {
    stop("basis_logE: base point must a square matrix.")
  }
  
  log_p = logm(p)
  idx = which(upper.tri(matrix(1, m, m), diag = TRUE), arr.ind = TRUE)
  num_basis = m * (m + 1) / 2
  res = array(NA, dim = c(num_basis, m, m))
  for (i in 1:num_basis) {
    E_i = matrix(0, m, m)
    row_idx = idx[i, 1]
    col_idx = idx[i, 2]
    
    if (row_idx == col_idx) {
      E_i[row_idx, col_idx] = 1.0
    } else {
      E_i[row_idx, col_idx] = 1.0 / sqrt(2)
      E_i[col_idx, row_idx] = 1.0 / sqrt(2)
    }
    
    res[i,,] = diff_explog(P = log_p, Q = E_i, type = "exp")
  }
  
  return (res)
}

#' @rdname basis
#' @export
basis.manifold_logE = function (mfd, p, ...) {
  basis_logE(p)
}

#' Evaluate the Riemannian metric at p (Log-Euclidean)
#' 
#' @param p an \eqn{m \times m} SPD matrix (base point)
#' @param v an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric matrices
#'          (tangent vectors at p)
#' @param w an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric matrices
#'          (tangent vectors at p)
#' 
#' @examples 
#' X = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' B = basis_logE(X)
#' Riem_metric_logE(X, B, B[1,,])
#' 
#' @export
Riem_metric_logE = function (p, v, w) {
  
  dim_v = dim(v)
  dim_w = dim(w)
  
  # Format v
  if (length(dim_v) == 2) {
    n_v = 1
    v = array(v, dim = c(1, dim_v))
  } else if (length(dim_v) == 3) {
    n_v = dim_v[1]
  } else {
    stop("Riem_metric_logE: dimension of v must be 2 or 3.")
  }
  
  # Format w
  if (length(dim_w) == 2) {
    n_w = 1
    w = array(w, dim = c(1, dim_w))
  } else if (length(dim_w) == 3) {
    n_w = dim_w[1]
  } else {
    stop("Riem_metric_logE: dimension of w must be 2 or 3.")
  }
  
  if (n_v != n_w) {
    if (n_v == 1) {
      v = array(rep(v[1,,], n_w), dim = c(dim_v, n_w))
      v = aperm(v, c(3, 1, 2))
      n = n_w
    } else if (n_w == 1) {
      w = array(rep(w[1,,], n_v), dim = c(dim_w, n_v))
      w = aperm(w, c(3, 1, 2))
      n = n_v
    } else {
      stop("Riem_metric_logE: number of tangent vectors in v and w must match or be 1.")
    }
  } else {
    n = n_v
  }
  
  m = dim(p)[1]
  if (length(dim(p)) != 2 || dim(p)[2] != m || dim(v)[2] != m || dim(v)[3] != m) {
    stop("Riem_metric_logE: input matrix dimensions must match (m x m).")
  }
  
  # --- Computation ---
  res = numeric(n)
  
  for (i in 1:n) {
    v_flat = diff_explog(P = p, Q = v[i,,], type = "log")
    w_flat = diff_explog(P = p, Q = w[i,,], type = "log")
    
    res[i] = sum(v_flat * w_flat)
  }
  
  if (n == 1) {
    return (res[1])
  }
  return(res)
}

#' @rdname Riem_metric
#' @export
Riem_metric.manifold_logE = function (mfd, p, v, w, ...) {
  Riem_metric_logE(p, v, w)
}

#' Compute the Riemannian Hessian vector action \eqn{H[v]} on the Log-Euclidean geometry
#' 
#' Evaluates the action of the Riemannian Hessian of \eqn{f(x) = 0.5 * d^2(x, \mu)}
#' on one or more tangent vectors v
#' 
#' @param x  base point where the Hessian is evaluated
#' @param p  target reference point
#' @param V  tangent vector(s) at x (\eqn{m \times m} or \eqn{n \times m \times m} arrays)
#' 
#' @export
Hess_logE = function (x, p, V) {
  return (V) # Hessian is identity
}

#' @rdname Hessian
#' @export
Hessian.manifold_logE = function (mfd, p, x, V, ...) {
  Hess_logE(x = x, p = p, V = V)
}

#' Fast Fréchet mean computation for log-Euclidean space
#' 
#' @param x an \eqn{n \times m \times m} array of SPD matrices
#' @param method either "specialized" (default), which uses specialized algorithm
#'               for log-Euclidean geometry, or "SGD" which uses generic Riemannian
#'               SGD
#' @return the Fréchet mean, as an \eqm{m \times m} matrix 
#' 
#' @examples 
#' mfd = manifold_logE()
#' x = array(NA, dim = c(10, 3, 3))
#' for (i in 1:10) x[i,,] = crossprod(matrix(rnorm(9), ncol = 3)) + diag(0.1, 3)
#' frechet_mean(mfd, x)
#' frechet_mean(mfd, x, method = "SGD", max.iter = 50)
#' 
#' @references 
#'   Arsigny, V., Fillard, P., Pennec, X., & Ayache, N. (2007). 
#'   Geometric means in a novel vector space structure on symmetric positive-definite matrices. 
#'   \emph{SIAM Journal on Matrix Analysis and Applications}, 29(1), 328--347. 
#' 
#' @export
frechet_mean.manifold_logE = function (mfd, x, method = c("specialized", "SGD"), ...) {
  method = match.arg(method)
  
  if (method == "SGD") {
    return (frechet_mean.default(mfd, x, ...))
  }
  
  if (length(dim(x)) == 2) {
    return (x)
  } else if (length(dim(x)) == 3) {
    n = dim(x)[1]
    m = dim(x)[2]
    if (n == 1) return (x[1,,])
  } else {
    stop("frechet_mean: the dimension x must be either 2 or 3 for log-Euclidean 
         geometry")
  }
  Log_avg = matrix(0, nrow = m, ncol = m)
  for (i in 1:n) {
    Log_avg = Log_avg + logm(x[i,,]) / n
  }
  mu_hat = expm(Log_avg)
  
  return (mu_hat)
}

#' @rdname tangent_to_vec
#' @export
tangent_to_vec.manifold_logE = function (mfd, v, p, E, ...) {
  mfd_dim = dim(E)[1]
  if (is.matrix(v)) {
    n = 1
    v = array(v, dim = c(1, dim(v)))
  } else if (length(dim(v)) == 3) {
    n = dim(v)[1]
  } else {
    stop("tangent_to_vec: dimension of v must be either 2 or 3")
  }
  
  res = array(NA, dim = c(n, mfd_dim))
  for (i in 1:n) {
    res[i,] = Riem_metric_logE(p = p, v = v[i,,], w = E)
  }
  if (n == 1) {
    res = res[1,]
  }
  
  return (res)
}

#' @rdname vec_to_tangent
#' @export
vec_to_tangent.manifold_logE = function (mfd, v_coord, p, E, ...) {
  if (is.vector(v_coord)) {
    n_v = 1
    v_coord = matrix(v_coord, nrow = 1)
  } else if (is.matrix(v_coord)) {
    n_v = nrow(v_coord)
  } else {
    stop("vec_to_tangent: v_coord must be either vector or matrix")
  }
  m = dim(p)[1]
  
  E_mat = t(matrix(aperm(E, c(2, 3, 1)), nrow = m * m))
  res_mat = v_coord %*% E_mat
  res = array(res_mat, dim = c(n_v, m, m))
  if (n_v == 1) {
    res = res[1,,]
  }
  
  return (res)
}








