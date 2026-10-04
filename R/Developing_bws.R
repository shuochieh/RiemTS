library(maotai)
library(expm)
library(deSolve)

#' Fast lyapunov solver when coefficient matrix is symmetric
#' 
#' Solves \eqn{AX + XA^{\top} = Q}.
#' When A is a symmetric matrix, a faster routine is used
#' 
#' @param A
#' @param Q
#' 
#' @examples 
#' A = crossprod(matrix(rnorm(500 * 500, sd = 0.5), ncol = 500)) + diag(0.5, 500)
#' X = matrix(rnorm(500 * 500, sd = 0.5), ncol = 500)
#' X = 0.5 * (X + t(X))
#' Q = A %*% X + X %*% A
#' 
#' norm(fast_lyapunov(A, Q) - X, "F")
#' 
#' @export
fast_lyapunov = function (A, Q) {
  if (!isSymmetric.matrix(A)) {
    return (lyapunov(A, Q))
  }
  
  temp = eigen(A)
  P = temp$vectors
  D = temp$values
  
  aux = outer(D, D, FUN = "+")
  M = t(P) %*% Q %*% P
  M = M / aux
  
  res = P %*% M %*% t(P)
  
  return (res)
}


#' Geodesic distance between two SPD matrices in Bures--Wasserstein
#' 
#' @noRd
geod_BWS_core = function (x, y) {
  # drop unnecessary dimension
  if (length(dim(x)) == 3 && dim(x)[1] == 1) {
    x = x[1,,]
  } 
  if (length(dim(y)) == 3 && dim(y)[1] == 1) {
    y = y[1,,]
  }
  
  temp = sqrtm(x)
  M = temp %*% y %*% temp
  M = 0.5 * (M + t(M)) # avoid asymmetry caused by numerical instability
  M = sqrtm(M)
  res = sum(diag(x)) + sum(diag(y)) - 2 * sum(diag(M))
  
  res = max(0, res) # avoid tiny negative values
  
  return (sqrt(res))
}

#' Geodesic distance between two (arrays of) points on the Bures--Wasserstein space
#' 
#' @param x an \eqn{m \times m} or \eqn{n \times m \times m} array of SPD matrix
#' @param y an \eqn{m \times m} or \eqn{n \times m \times m} array of SPD matrix
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' y = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' geod_BWS(x, y)
#' 
#' @export
geod_BWS = function (x, y) {
  if (length(dim(x)) == 3) {
    nx = dim(x)[1]
  } else if (length(dim(x)) == 2) {
    nx = 1
  } else {
    stop("geod_BWS: input array dimension must be 2 or 3.")
  }
  mx = dim(x)[2]
  
  if (length(dim(y)) == 3) {
    ny = dim(y)[1]
  } else if (length(dim(y)) == 2) {
    ny = 1
  } else {
    stop("geod_BWS: input array dimension must be 2 or 3.")
  }
  my = dim(y)[2]
  
  if (nx != 1 && ny != 1 && nx != ny) {
    stop("geod_BWS: number of input matrices must match or equal to 1; otherwise broadcasting is undefined.")
  } 
  if (mx != my) {
    stop("geod_BWS: input matrix dimensions must match.")
  }
  
  m = mx
  
  if (nx == 1 && ny == 1) {
    return (geod_BWS_core(x, y))
  } else if (nx == 1) {
    res = rep(NA, ny)
    for (i in 1:ny) {
      res[i] = geod_BWS_core(x, y[i,,])
    }
  } else if (ny == 1) {
    res = rep(NA, nx)
    for (i in 1:nx) {
      res[i] = geod_BWS_core(x[i,,], y)
    }
  } else {
    res = rep(NA, nx)
    for (i in 1:nx) {
      res[i] = geod_BWS_core(x[i,,], y[i,,])
    }
  }
  
  return (res)
}

#' @rdname geod
#' @export
geod.manifold_BWS = function (mfd, x, y, ...) {
  geod_BWS(x, y)
}

#' Exponential map on the Bures--Wasserstein geometry
#' 
#' @noRd
Exp_BWS_core = function (z, x) {
  if (length(dim(x)) == 3 && dim(x)[1] == 1) {
    x = x[1,,]
  }
  
  L = lyapunov_fast(x, z)
  d = dim(z)[1]
  L = L + diag(1, d)
  res = L %*% x %*% L
  
  return (res)
}

#' Exponential map on the Bures--Wasserstein geometry
#' 
#' @param z an \eqn{m \times m}  or \eqn{n \times m \times m} array of tangent vectors (identified as symmetric matrices)
#' @param x an \eqn{m \times m} SPD matrix (as the base point)
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)
#' z = crossprod(matrix(rnorm(9), 3, 3)) * 0.1
#' Exp_BWS(z, x)
#'
#' @export
Exp_BWS = function (z, x) {
  if (length(dim(z)) == 3) {
    n = dim(z)[1]
    was_matrix = FALSE
  } else if (length(dim(z)) == 2) {
    n = 1
    z = array(z, dim = c(1, dim(z)))
    was_matrix = TRUE
  } else {
    stop("Exp_BWS: number of dimensions of z must be either 2 or 3.")
  }
  
  m = dim(z)[2]
  if (dim(z)[3] !=m || !all(dim(x) == c(m, m))) {
    stop("Exp_BWS: input matrix dimension must match (m x m).")
  }
  
  res = array(NA, dim = c(n, m, m))
  for (i in 1:n) {
    res[i,,] = Exp_BWS_core(z = z[i,,], x = x)
  }
  if (was_matrix) {
    res = res[1,,]
  }
  
  return (res)
}

#' @export
Exp_mfd.manifold_BWS = function (mfd, p, v, ...) {
  Exp_BWS(z = v, x = p)
}

#' Logarithm map on the Bures--Wasserstein geometry
#' 
#' @noRd
Log_BWS_core = function (x, y, check_cut = FALSE) {
  if (length(dim(x)) == 3 && dim(x)[1] == 1) {
    x = x[1,,]
  }
  
  if (check_cut) {
    # check if y is in the injectivity radius
    test_y = fast_lyapunov(x, y)
    test_y = test_y + diag(1, nrow(y))
    eig = eigen(test_y, symmetric = T, only.values = T)$values
    if (!all(eig > 0)) {
      stop("Log_BWS_core: y outside injectivity radius")
    }
  }
  
  res = sqrtm(y %*% x)
  res = res + t(res) - 2 * x
  
  return (res)
}

#' Logarithm map on the Bures--Wasserstein geometry
#' 
#' @param x an \eqn{m \times m} SPD matrix (as the base point)
#' @param y an \eqn{m \times m}  or \eqn{n \times m \times m} array of SPD matrices
#' @param check_cut whether to perform checks for injectivity radius
#' 
#' @examples 
#' x = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' y = crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' Log_BWS(x, y)
#'  
#' @export
Log_BWS = function (x, y, check_cut = FALSE) {
  if (length(dim(y)) == 3) {
    n = dim(y)[1]
    was_matrix = FALSE
  } else if (length(dim(y)) == 2) {
    n = 1
    y = array(y, dim = c(1, dim(y)))
    was_matrix = TRUE
  } else {
    stop("Log_BWS: number of dimensions of y must be either 2 or 3.")
  }
  
  m = dim(y)[2]
  if (dim(y)[3] != m || !all(dim(x) == c(m, m))) {
    stop("Log_BWS: input matrix dimension must match (m x m).")
  }
  
  res = array(NA, dim = c(n, m, m))
  for (i in 1:n) {
    res[i,,] = Log_BWS_core(x = x, y = y[i,,])
  }
  if (was_matrix) {
    res = res[1,,]
  }
  
  return (res)
}

#' @rdname Log_mfd
#' @export
Log_mfd.manifold_BWS = function (mfd, p, q, ...) {
  Log_BWS(x = p, y = q, ...)
}

#' @noRd
Christoffel_BWS_core = function (Sigma, X, Y) {
  Lx = lyapunov(Sigma, X)
  Ly = lyapunov(Sigma, Y)
  
  res = (Sigma %*% Ly %*% Lx) + (Ly %*% Lx %*% Sigma) - (Lx %*% Y) - (Ly %*% X)
  res = 0.5 * (res + t(res))
  
  return (res)
}

#' @noRd
Christoffel_BWS = function (t, U, param) {
  
  Sigma1 = param$Sigma1
  Sigma2 = param$Sigma2
  p = dim(Sigma1)[1]
  U = matrix(U, p, p)
  
  A = Sigma1 %*% Sigma2
  B = Sigma2 %*% Sigma1
  A = sqrtm(A)
  B = sqrtm(B)
  
  Sigma_t = (1 - t)^2 * Sigma1 + t^2 * Sigma2 + t * (1 - t) * (A + B)
  Sigma_dot = -2 * Sigma1 + A + B + 2 * t * (Sigma1 + Sigma2 - A - B)
  
  res = -as.vector(Christoffel_BWS_core(Sigma_t, Sigma_dot, U))
  return (list(res))
}

#' Computes parallel transport in the Bures-Wasserstein space along geodesic
#' 
#' @param Sigma1 starting point
#' @param Sigma2 end point
#' @param V tangent vector, identified as a symmetric matrix, at the starting point
#' 
#' @noRd
pt_bws_core = function (Sigma1, Sigma2, V, method = "adams") {
  p = dim(Sigma1)[1]
  times = seq(0, 1, length.out = 101)
  sol = ode(y = as.vector(V),
            times = times,
            func = Christoffel_BWS,
            parms = list("Sigma1" = Sigma1, "Sigma2" = Sigma2),
            method = method)
  
  res = matrix(sol[101, -1], p, p)
  return (res)
}

#' Parallel transport along geodesic on the Bures--Wasserstein geometry
#' 
#' @param p an \eqn{m \times m} SPD matrix (start)
#' @param q an \eqn{m \times m} SPD matrix (end)
#' @param x an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric
#'          matrices, identified as tangent vectors at p
#'          
#' @examples 
#' p = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' q = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' v = matrix(rnorm(9), ncol = 3)
#' v = 0.5 * (v + t(v))
#' pt_BWS(p, q, v)
#' 
#' @export
pt_BWS = function (p, q, x, method = "adams") {
  if (length(dim(x)) == 3) {
    n = dim(x)[1]
    was_matrix = FALSE
  } else if (length(dim(x)) == 2) {
    n = 1
    x = array(x, dim = c(1, dim(x)))
    was_matrix = TRUE
  } else {
    stop("pt_BWS: number of dimensions of x must be either 2 or 3.")
  }
  
  m = dim(x)[2]
  res = array(NA, dim = c(n, m, m))
  for (i in 1:n) {
    res[i,,] = pt_bws_core(Sigma1 = p, Sigma2 = q, V = x[i,,], method = method)
  }
  if (was_matrix) {
    res = res[1,,]
  }
  
  return (res)
}

#' @rdname ptransport
#' @export
ptransport.manifold_BWS = function (mfd, from, to, v, ...) {
  pt_BWS(p = from, q = to, x = V, ...)
}

#' A basis for the tangent space at p (SPD matrix in the Bures--Wasserstein geometry)
#'
#' @param p an \eqn{m \times m} SPD matrix (base point)
#' 
#' @examples 
#' p = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' basis_BWS(p)
#' 
#' @export
basis_BWS = function (p) {
  d = dim(p)[1]
  
  model = eigen(p)
  lambdas = model$values
  P = model$vectors
  
  E = array(NA, dim = c(d * (d + 1) / 2, d, d))
  # E_lyapunov = array(NA, dim = c(d * (d + 1) / 2, d, d))
  counter = 0
  for (i in 1:d) {
    for (j in i:d) {
      counter = counter + 1
      S = matrix(0, ncol = d, nrow = d)
      # S_tilde = matrix(0, ncol = d, nrow = d)
      if (i == j) {
        S[i,j] = sqrt(2 * (lambdas[i] + lambdas[j]))
        # S_tilde[i,j] = 1 / sqrt(lambdas[i])
      } else {
        S[i,j] = sqrt(lambdas[i] + lambdas[j])
        S[j,i] = sqrt(lambdas[i] + lambdas[j])
        
        # S_tilde[i,j] = 1 / sqrt(lambdas[i] + lambdas[j])
        # S_tilde[j,i] = 1 / sqrt(lambdas[i] + lambdas[j])
      }
      E[counter,,] = P %*% S %*% t(P)
      # E_lyapunov[counter,,] = P %*% S_tilde %*% t(P)
    }
  }
  
  # return (list("E" = E, "E_lyapunov" = E_lyapunov))
  return (E)
}

#' @rdname basis
#' @export
basis.manifold_BWS = function (mfd, p, ...) {
  basis_BWS(p)
}

#' Evaluate the Riemannian metric at p (Bures--Wasserstein)
#' 
#' @param p an \eqn{m \times m} SPD matrix (base point)
#' @param v an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric matrices
#'          (tangent vectors at p)
#' @param w an \eqn{m \times m} or \eqn{n \times m \times m} array of symmetric matrices
#'          (tangent vectors at p)
#' 
#' @examples 
#' X = crossprod(matrix(rnorm(9), ncol = 3)) + diag(1, 3)
#' B = basis_BWS(X)
#' Riem_metric_BWS(X, B, B[1,,])
#' 
#' @export
Riem_metric_BWS = function (p, v, w) {
  
  dim_v = dim(v)
  dim_w = dim(w)
  
  # Format v
  if (length(dim_v) == 2) {
    n_v = 1
    v = array(v, dim = c(1, dim_v))
  } else if (length(dim_v) == 3) {
    n_v = dim_v[1]
  } else {
    stop("Riem_metric_BWS: dimension of v must be 2 or 3.")
  }
  
  # Format w
  if (length(dim_w) == 2) {
    n_w = 1
    w = array(w, dim = c(1, dim_w))
  } else if (length(dim_w) == 3) {
    n_w = dim_w[1]
  } else {
    stop("Riem_metric_BWS: dimension of w must be 2 or 3.")
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
      stop("Riem_metric_BWS: number of tangent vectors in v and w must match or be 1.")
    }
  } else {
    n = n_v
  }
  
  m = dim(p)[1]
  if (length(dim(p)) != 2 || dim(p)[2] != m || dim(v)[2] != m || dim(v)[3] != m) {
    stop("Riem_metric_BWS: input matrix dimensions must match (m x m).")
  }
  
  # --- Computation ---
  if (n_v == 1) {
    v_lyapunov = fast_lyapunov(p, v[1,,])
    w_mat = matrix(aperm(w, c(2, 3, 1)), nrow = m * m)
    res = 0.5 * c(as.vector(v_lyapunov) %*% w_mat)
    
    return (res)
  } else if (n_w == 1) {
    w_lyapunov = fast_lyapunov(p, w[1,,])
    v_mat = matrix(aperm(v, c(2, 3, 1)), nrow = m * m)
    res = 0.5 * c(as.vector(w_lyapunov) %*% v_mat)
    
    return (res)
  }

  res = numeric(n)
  for (i in 1:n) {
    v_lyapunov = lyapunov(p, v[i,,])
    res[i] = 0.5 * sum(diag(v_lyapunov %*% w[i,,]))
    # if (n_v == 1) {
    #   res[i] = 0.5 * sum(diag(v_lyapunov %*% w[i,,]))
    # } else if (n_w == 1) {
    #   res[i] = 0.5 * sum(diag(w_lyapunov %*% v[i,,]))
    # } else {
    #   v_lyapunov = lyapunov(p, v[i,,])
    #   res[i] = 0.5 * sum(diag(v_lyapunov %*% w[i,,]))
    # }
  }
  
  if (n == 1) {
    return (res[1])
  }
  return (res)
}

#' @rdname Riem_metric
#' @export
Riem_metric.manifold_BWS = function (mdf, p, v, w, ...) {
  Riem_metric_BWS(p = p, v = v, w = w)
}

#' Compute the Riemannian Hessian vector action \eqn{H[v]} on the Bures--Wasserstein geometry
#' 
#' Evaluates the action of the Riemannian Hessian of \eqn{f(x) = 0.5 * d^2(x, \mu)}
#' on one or more tangent vectors v
#' 
#' @param x  base point where the Hessian is evaluated
#' @param p  target reference point
#' @param V  tangent vector(s) at x (\eqn{m \times m} or \eqn{n \times m \times m} arrays)
#' 
# To be done...


#' Fast Frèchet mean computation for Bures--Wasserstein space
#' 
#' @param x an \eqn{n \times m \times m} array of SPD matrices
#' @param method either "specialized" (default), which uses specialized algorithm
#'               (if supported), or "SGD" which uses generic Riemannian SGD
#' @param max.iter maximum number of iterations
#' @param tol convergence tolerance
#' @param verbose whether to print loss along iterations
#' 
#' @return the Fréchet mean, as an \eqm{m \times m} matrix 
#' 
#' @examples 
#' mfd = manifold_bws()
#' x = array(NA, dim = c(10, 3, 3))
#' for (i in 1:10) x[i,,] = crossprod(matrix(rnorm(9), ncol = 3)) + diag(0.1, 3)
#' frechet_mean(mfd, x, verbose = TRUE)
#' frechet_mean(mfd, x, method = "SGD", verbose = TRUE)
#' 
#' @references 
#'  Álvarez-Esteban, P. C., del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2016). 
#'  A fixed-point approach to barycenters in Wasserstein space.
#'  \emph{Journal of Mathematical Analysis and Applications}, 441(2), 744--762.
#' 
#' @export
frechet_mean.manifold_BWS = function (mfd, x, method = c("specialized", "SGD"), 
                                      max.iter = 100, tol = 1e-4, verbose = FALSE,
                                      ...) {
  method = match.arg(method)
  
  if (method == "SGD") {
    return (frechet_mean.default(mfd, x, max.iter = max.iter, tol = tol, 
                                 verbose = verbose, ...))
  }
  
  if (length(dim(x)) == 2) {
    return (x)
  } else if (length(dim(x)) == 3) {
    n = dim(x)[1]
    m = dim(x)[2]
    if (n == 1) return (x[1,,])
  } else {
    stop("frechet_mean: the dimension of x must be either 2 or 3 for Bures--Wasserstein geometry")
  }
  
  S = x[sample(n, 1),,]
  for (zz in 1:max.iter) {
    S_eigen = eigen(S, symmetric = TRUE)
    S_half = S_eigen$vectors %*% sqrt(diag(pmax(S_eigen$values, 0))) %*% t(S_eigen$vectors)
    S_ihalf = S_eigen$vectors %*% diag(1 / sqrt(pmax(S_eigen$values, 0))) %*% t(S_eigen$vectors)
    
    M = matrix(0, nrow = m, ncol = m)
    for (i in 1:n) {
      A = S_half %*% x[i,,] %*% S_half
      A_eigen = eigen(A, symmetric = TRUE)
      
      A_half = A_eigen$vectors %*% sqrt(diag(pmax(A_eigen$values, 0))) %*% t(A_eigen$vectors)
      
      M = M + (A_half / n)
    }
    
    S_new = S_ihalf %*% (M %*% M) %*% S_ihalf
    
    
    loss = mean(sum((S_new - S)^2)) # normalized Frobenius 
    if (verbose) {
      cat("frechet_mean:", class(mfd)[1], "iter", zz, "(fixed point) loss:",
          sprintf("%.4f", loss), "\n")
    }
    
    if (loss < tol) {
      S = S_new
      break
    }
    S = S_new
  }
  
  return (S)
}

#' @rdname tangent_to_vec 
#' @export
tangent_to_vec.manifold_BWS = function (mfd, v, p, E, ...) {
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
    res[i,] = Riem_metric_BWS(p = p, v = v[i,,], w = E)
  }
  if (n == 1) {
    res = res[1,]
  }
  
  return (res)
}


#' @rdname vec_to_tangent
#' @export
vec_to_tangent.manifold_BWS = function (mfd, v_coord, p, E, ...) {
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




