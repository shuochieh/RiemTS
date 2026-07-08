################################################################################
#### Bures-Wasserstein manifold: core math and S3 methods #####################
################################################################################

#' Geodesic distance between two SPD matrices in Bures-Wasserstein
#' @export
geod_BWS_core = function(X, Y) {
  temp = sqrtm(X)
  temp = sqrtm(temp %*% Y %*% temp)
  res = sum(diag(X)) + sum(diag(Y)) - 2 * sum(diag(temp))

  return(sqrt(res))
}

#' @export
geod_core.manifold_bws = function(mfd, x, y, ...) geod_BWS_core(x, y)

#' Exponential map on the Bures-Wasserstein manifold
#' @export
Exp_BWS_core = function(X, M) {
  L = lyapunov(M, X)
  diag(L) = diag(L) + 1
  output = L %*% M %*% L

  return(output)
}

#' @export
Exp_core.manifold_bws = function(mfd, v, mu, ...) Exp_BWS_core(v, mu)

#' Log map on the Bures-Wasserstein manifold
#' @export
Log_BWS_core = function(X, M) {
  # check if X is within the injectivity radius
  test_X = lyapunov(M, X)
  diag(test_X) = diag(test_X) + 1
  eig = eigen(test_X, symmetric = TRUE, only.values = TRUE)$values
  if (!all(eig > 0)) {
    stop("Log_BWS: X outside injectivity radius")
  }

  output = sqrtm(X %*% M)
  output = output + sqrtm(M %*% X) - 2 * M

  return(output)
}

#' @export
Log_core.manifold_bws = function(mfd, x, mu, ...) Log_BWS_core(x, mu)

#' Computes a canonical basis for the tangent space at a point in BWS
#'
#' @param Sigma Base point
#' @return An array of orthonormal basis (E), and an array of E passed through the Lyapunov operator
#' @export
tan_basis_bws = function(Sigma) {
  p = dim(Sigma)[1]

  model = eigen(Sigma)
  lambdas = model$values
  P = model$vectors

  E = array(NA, dim = c(p * (p + 1) / 2, p, p))
  E_lyapunov = array(NA, dim = c(p * (p + 1) / 2, p, p))
  counter = 0
  for (i in 1:p) {
    for (j in i:p) {
      counter = counter + 1
      S = matrix(0, ncol = p, nrow = p)
      S_tilde = matrix(0, ncol = p, nrow = p)
      if (i == j) {
        S[i, j] = sqrt(2 * (lambdas[i] + lambdas[j]))
        S_tilde[i, j] = 1 / sqrt(lambdas[i])
      } else {
        S[i, j] = sqrt(lambdas[i] + lambdas[j])
        S[j, i] = sqrt(lambdas[i] + lambdas[j])

        S_tilde[i, j] = 1 / sqrt(lambdas[i] + lambdas[j])
        S_tilde[j, i] = 1 / sqrt(lambdas[i] + lambdas[j])
      }
      E[counter, , ] = P %*% S %*% t(P)
      E_lyapunov[counter, , ] = P %*% S_tilde %*% t(P) # L_{mu_hat}(E)
    }
  }

  return(list("E" = E, "E_lyapunov" = E_lyapunov))
}

#' Express the log-mapped data in E coordinates
#'
#' Express the tangent vector Log_M(x) in the canonical orthonormal basis
#' returned by tan_basis_bws(). Thin wrapper around Log_map()/coords_from_basis().
#'
#' @export
log_vec_construct = function(x, M, E_lyapunov = NULL) {
  if (is.null(E_lyapunov)) {
    E_lyapunov = tan_basis_bws(M)$E_lyapunov
  }

  log_x = Log_map(manifold_bws(), x, M)
  return(coords_from_basis(E_lyapunov, log_x))
}

Christoffel_BWS_core = function(Sigma, X, Y) {
  Lx = lyapunov(Sigma, X)
  Ly = lyapunov(Sigma, Y)

  res = (Sigma %*% Ly %*% Lx) + (Ly %*% Lx %*% Sigma) - (Lx %*% Y) - (Ly %*% X)
  res = 0.5 * (res + t(res))

  return(res)
}

Christoffel_BWS = function(t, U, param) {
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
  return(list(res))
}

#' Parallel transport in the Bures-Wasserstein space along the geodesic
#'
#' @param Sigma1 starting point
#' @param Sigma2 end point
#' @param V tangent vector, identified as a symmetric matrix, at the starting point
pt_bws_core = function(Sigma1, Sigma2, V, method = "adams") {
  p = dim(Sigma1)[1]
  times = seq(0, 1, length.out = 101)
  sol = ode(
    y = as.vector(V),
    times = times,
    func = Christoffel_BWS,
    parms = list("Sigma1" = Sigma1, "Sigma2" = Sigma2),
    method = method
  )

  res = matrix(sol[101, -1], p, p)
  return(res)
}

#' @export
parallel_transport.manifold_bws = function(mfd, from, to, v, method = "adams", ...) {
  pt_bws_core(from, to, v, method)
}
