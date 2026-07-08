################################################################################
#### Utility functions for analysis of data in the Grassmannian ################
################################################################################
library(Matrix)

#' computes the geodesic distance between two data point in Gr(d, p)
#' 
geod_gr_core = function (x, y, p = NULL) {
  if (!any(dim(x) == dim(y))) {
    stop("goed_gr_core:x, y has different ambient dimensions")
  }
  
  if (is.null(p)) {
    p = round(sum(eigen(x, symmetric = TRUE, only.values = TRUE)$values > 0.5))
    # assumed to be equal for x and y
  }
  
  lambda = eigen(x %*% y, symmetric = TRUE, only.values = TRUE)$values
  lambda = Re(lambda)
  lambda = sort(lambda, decreasing = TRUE)[1:p]
  lambda = pmin(1, pmax(lambda, 0))
  
  theta = acos(sqrt(lambda))
  
  return (sqrt(sum(theta^2)))
}

#' computes the geodesic distance between two data point(s) in Gr(d, p)
#' 
#' @param x a projection matrix (d by d) or an (n by d by d) array of projectors
#' @param y same as x
#' @param p rank of the projectors in the Grassmannian
#' 
geod_gr = function (x, y, p = NULL) {
  
  if ((length(dim(x)) == 2) && (length(dim(y)) == 2)) {
    
    res = geod_gr_core(x, y, p = p)
    
  } else if ((length(dim(x)) == 2) && (length(dim(y)) == 3)) {
    
    res = rep(NA, dim(y)[1])
    for (i in 1:dim(y)[1]) {
      res[i] = geod_gr_core(x, y[i,,], p = p)
    }
    
  } else if ((length(dim(x)) == 3) && (length(dim(y)) == 2)) {
      
    res = rep(NA, dim(x)[1])
    for (i in 1:dim(x)[1]) {
      res[i] = geod_gr_core(x[i,,], y, p = p)
    }
    
  } else if ((length(dim(x)) == 3) && (length(dim(y)) == 3)) {
    
    res = rep(NA, dim(x)[1])
    if (dim(x)[1] != dim(y)[1]) {
      stop("geod_gr: number of projectors in x and y differs")
    }
    for (i in 1:dim(x)[1]) {
      res[i] = geod_gr_core(x[i,,], y[i,,], p = p)
    }
  }
  
  return (res)
}

#' exponential map on Gr(d, p)
#' 
#' @param D A projector-perspective tangent element (a symmetric matrix)
#' @param P A projector (base of the exponential)
#' 
Exp_gr_core = function (D, P) {
  Comm = D %*% P - P %*% D
  
  exp_C = expm(Comm)
  exp_mC = t(exp_C)
  # exp_mC = solve(exp_C)
  
  temp = exp_C %*% P %*% exp_mC
  res = (temp + t(temp)) / 2 # numerically ensure it is symmetric
  
  return (as.matrix(res))
}

#' exponential map on Gr(d, p)
#' 
#' @param D A projector-perspective tangent element (a symmetric matrix) or array of such elements
#' @param P A projector (base of the exponential)
#' 
Exp_gr = function (D, P) {
  if (length(dim(D)) == 2) {
    res = Exp_gr_core(D, P)
  } else if (length(dim(D)) == 3) {
    res = array(NA, dim = dim(D))
    if (!any(dim(D)[c(2,3)] == dim(P))) {
      stop ("Exp_gr: dimension of D and P not aligned")
    }
    for (i in 1:dim(D)[1]) {
      res[i,,] = Exp_gr_core(D[i,,], P)
    }
  } else {
    stop("Exp_gr: dimension of D not supported")
  }
  
  return (res)
}

#' Riemannian log map on Gr(d, p)
#' 
#' @param P base
#' @param x argument
#' @param p rank of the projectors in the Grassmannian
#' 
Log_gr_core = function (P, x, p = NULL) {
  if (is.null(p)) {
    p = round(sum(eigen(P, symmetric = TRUE, only.values = TRUE)$values > 0.5))
  }
  
  U = eigen(P, symmetric = TRUE)$vectors[,1:p]
  Y = eigen(x, symmetric = TRUE)$vectors[,1:p]
  
  temp = svd(t(Y) %*% U)
  Q_tilde = temp$u
  S_tilde = temp$d
  R_tilde = temp$v
  
  R = R_tilde[,rev(1:ncol(R_tilde))]
  S = S_tilde[rev(1:length(S_tilde))]
  S = pmin(S, 1)
  S_hat = sqrt((1 - S^2))
  
  if (sum(S_hat) < 1e-8) {
    return (x)
  }
  
  S_hat_dup = S_hat
  S_hat_dup[which(abs(S_hat) < 1e-8)] = 1 # for numerical stability (these entries will be killed by asin(S_hat) later)
  
  Q = Q_tilde[,rev(1:ncol(Q_tilde))]
  Q_hat = (diag(1, nrow(P)) - P) %*% Y %*% Q %*% diag(S_hat_dup^(-1), nrow = p)
  
  Sigma = asin(S_hat)
  D = Q_hat %*% diag(Sigma, nrow = p) %*% t(R)
  
  return (U %*% t(D) + D %*% t(U))
}

#' Riemannian log map on Gr(d, p)
#' 
#' @param P base (a fixed projector)
#' @param x argument (an array of projectors)
#' @param p rank of the projectors in the Grassmannian
#' 
Log_gr = function (P, x, p = NULL) {
  if (length(dim(x)) == 2) {
    res = Log_gr_core(P, x, p = p)
  } else if (length(dim(x)) == 3) {
    n = dim(x)[1]
    res = array(NA, dim = dim(x))
    for (i in 1:n) {
      res[i,,] = Log_gr_core(P, x[i,,], p = p)
    }
  } else {
    stop("Log_gr: x must 3-d array or a matrix")
  }
  
  return (res)
}

#' Generate a set of ONB at a given P
#' 
#' @param P base
#' @param p rank of the projectors in the Grassmannian
#' 
basis_gr = function (P, p = NULL) {
  
  d = nrow(P)
  temp = eigen(P, symmetric = TRUE)
  
  if (is.null(p)) {
    p = round(sum(temp$values > 0.5))
  }
  
  U = temp$vectors[,1:p,drop = FALSE]
  U_perp = temp$vectors[,(p + 1):d,drop = FALSE]
  
  m = (d - p) * p  # intrinsic dimension
  
  basis_array = array(NA, dim = c(m, d, d))
  idx = 1
  for (a in 1:(d - p)) {
    for (b in 1:p) {
      E = matrix(0, d - p, p)
      E[a, b] = 1
      Xi = U_perp %*% E %*% t(U) + U %*% t(E) %*% t(U_perp)
      
      basis_array[idx,,] = Xi
      
      idx = idx + 1
    }
  }
  
  return (basis_array)
}

#' get coordinates with respect to certain basis
#' 
#' @param basis 
#' @param x (d by d) or (n by d by d)
#' 
coord_gr = function (basis, x) {
  
  m = dim(basis)[1]
  
  if (length(dim(x)) == 3) {
    n = dim(x)[1]
    res = array(NA, dim = c(n, m))
  } else if (length(dim(x)) == 2) {
    n = 1
    res = rep(NA, m)
  } else {
    stop("coord_gr: unsupported x dimension")
  }
  
  for (i in 1:m) {
    if (length(dim(x)) == 2) {
      res[i] = 0.5 * t(c(basis[i,,])) %*% c(x)
    } else {
      for (j in 1:n) {
        res[j,i] = 0.5 * t(c(basis[i,,])) %*% c(x[j,,])
      }
    }
  }
  
  return (res)
}

#' get tangent vector from coordinates (in projector perspectives)
#' 
coord2tan = function (basis, coord) {
  m = dim(basis)[1]
  d = dim(basis)[2]
  res = matrix(0, nrow = d, ncol = d)
  for (i in 1:m) {
    res = res + coord[i] * basis[i,,]
  }
  
  return (res)
}

#' computes the Frechet mean on the Grassmannian
#' 
#' @param x (n by d by d) array of data
#' @param p rank of the projectors in the Grassmannian
#' 
mean_on_gr = function (x, p, tau = 0.1, tol = 1e-8, max.iter = 1000,
                       batch_por = 1.0, init = NULL, verbose = FALSE) {
  n = dim(x)[1]
  if (n == 1) {
    return (x[1,,])
  }
  
  if (is.null(init)) {
    mu = x[sample(n, 1),,]
  } else {
    mu = init
  }
  
  for (i in 1:max.iter) {
    basis = basis_gr(mu, p = p)
    
    x_batch = x[sample(n, floor(batch_por * n)),,]
    temp = Log_gr(mu, x_batch, p = p)
    coord = coord_gr(basis, temp)
    grad_coord = colMeans(coord)
    grad = coord2tan(basis, grad_coord)
    
    mu_new = Exp_gr(tau * grad, mu)
    
    loss = mean(geod_gr(x, mu_new))
    if (verbose) {
      cat("iter", i, ":", loss, "\n")
    }
    
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
  }
  
  return (mu)
}

#' parallel transport via geodesics on the Grassmannian
#' 
#' @param P starting point
#' @param W tangent vector associated with the end point (Log_start End)
#' @param V tangent vector at P1 (a matrix)
#' 
pt_gr_core = function (P, W, V) {
  Comm = W %*% P - P %*% W
  E = expm(Comm)
  # E_m = solve(E)
  E_m = t(E)
  
  return (as.matrix(E %*% V %*% E_m))
}




