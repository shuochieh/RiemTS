################################################################################
#### Grassmannian manifold: core math and S3 methods ###########################
################################################################################

#' Geodesic distance between two projectors in Gr(d, p)
#' @export
geod_gr_core = function(x, y, p = NULL) {
  if (!any(dim(x) == dim(y))) {
    stop("geod_gr_core: x, y have different ambient dimensions")
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

  return(sqrt(sum(theta^2)))
}

#' @export
geod_core.manifold_grassmann = function(mfd, x, y, ...) geod_gr_core(x, y, p = mfd$p)

#' Exponential map on Gr(d, p)
#'
#' @param D A projector-perspective tangent element (a symmetric matrix)
#' @param P A projector (base of the exponential)
#' @export
Exp_gr_core = function(D, P) {
  Comm = D %*% P - P %*% D

  exp_C = expm(Comm)
  exp_mC = t(exp_C)

  temp = exp_C %*% P %*% exp_mC
  res = (temp + t(temp)) / 2 # numerically ensure it is symmetric

  return(as.matrix(res))
}

#' @export
Exp_core.manifold_grassmann = function(mfd, v, mu, ...) Exp_gr_core(v, mu)

#' Riemannian log map on Gr(d, p)
#'
#' @param P base
#' @param x argument
#' @param p rank of the projectors in the Grassmannian
#' @export
Log_gr_core = function(P, x, p = NULL) {
  if (is.null(p)) {
    p = round(sum(eigen(P, symmetric = TRUE, only.values = TRUE)$values > 0.5))
  }

  U = eigen(P, symmetric = TRUE)$vectors[, 1:p]
  Y = eigen(x, symmetric = TRUE)$vectors[, 1:p]

  temp = svd(t(Y) %*% U)
  Q_tilde = temp$u
  S_tilde = temp$d
  R_tilde = temp$v

  R = R_tilde[, rev(1:ncol(R_tilde))]
  S = S_tilde[rev(1:length(S_tilde))]
  S = pmin(S, 1)
  S_hat = sqrt((1 - S^2))

  if (sum(S_hat) < 1e-8) {
    return(x)
  }

  S_hat_dup = S_hat
  S_hat_dup[which(abs(S_hat) < 1e-8)] = 1 # numerical stability (killed by asin(S_hat) later)

  Q = Q_tilde[, rev(1:ncol(Q_tilde))]
  Q_hat = (diag(1, nrow(P)) - P) %*% Y %*% Q %*% diag(S_hat_dup^(-1), nrow = p)

  Sigma = asin(S_hat)
  D = Q_hat %*% diag(Sigma, nrow = p) %*% t(R)

  return(U %*% t(D) + D %*% t(U))
}

#' @export
Log_core.manifold_grassmann = function(mfd, x, mu, ...) Log_gr_core(mu, x, p = mfd$p)

#' Generate a set of ONB at a given P
#'
#' @param P base
#' @param p rank of the projectors in the Grassmannian
#' @export
basis_gr = function(P, p = NULL) {
  d = nrow(P)
  temp = eigen(P, symmetric = TRUE)

  if (is.null(p)) {
    p = round(sum(temp$values > 0.5))
  }

  U = temp$vectors[, 1:p, drop = FALSE]
  U_perp = temp$vectors[, (p + 1):d, drop = FALSE]

  m = (d - p) * p # intrinsic dimension

  basis_array = array(NA, dim = c(m, d, d))
  idx = 1
  for (a in 1:(d - p)) {
    for (b in 1:p) {
      E = matrix(0, d - p, p)
      E[a, b] = 1
      Xi = U_perp %*% E %*% t(U) + U %*% t(E) %*% t(U_perp)

      basis_array[idx, , ] = Xi

      idx = idx + 1
    }
  }

  return(basis_array)
}

#' Parallel transport via geodesics on the Grassmannian
#'
#' @param P starting point
#' @param W tangent vector associated with the end point (Log_start End)
#' @param V tangent vector at P (a matrix)
pt_gr_core = function(P, W, V) {
  Comm = W %*% P - P %*% W
  E = expm(Comm)
  E_m = t(E)

  return(as.matrix(E %*% V %*% E_m))
}

#' @export
parallel_transport.manifold_grassmann = function(mfd, from, to, v, ...) {
  W = Log_gr_core(from, to, p = mfd$p)
  pt_gr_core(from, W, v)
}
