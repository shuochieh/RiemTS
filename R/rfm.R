# Riemannian Factor Model
library(RSpectra)

#' Factor model of Lam, Yao, and Bathia (2011)
#' 
#' @param x data
#' @param r number of factors (default: NULL, in which case the ratio estimator
#'          is used to select it)
#' @param h number of lags to use in construction of the autocovariance-based 
#'          auxiliary matrix (default: 1)
#' @param demean whether to demean the data first (default: TRUE)
#' @param R maximum number of factors to consider in the ratio estimator (default:
#'          NULL, in which case p/2 is used, where p is the dimension of data)
#' @return A `factor_model` object containing
#' * `V` estimated loading matrix
#' * `f_hat` estimated factor process
#' * `e_hat` residuals
#' * `fitted.val` fitted values (pseudo-predictions)
#' * `r_hat` number of factors used
#' * `mean` means of each variable
#' 
#' @references 
#'   LAM, C., YAO, Q., & BATHIA, N. (2011). 
#'   Estimation of latent factors for high-dimensional time series. 
#'   Biometrika, 98(4), 901--918. 
#' 
#' @keywords internal
#' @export
LYB_fm = function (x, r = NULL, h = 1, demean = TRUE, R = NULL) {
  n = nrow(x)
  p = ncol(x)
  if (demean) {
    means = colMeans(x)
    x = t(t(x) - means)
  } else {
    means = rep(0, p)
  }
  
  # compute the auxiliary positive definite matrix
  L = 0
  H = h + 1
  for (i in 1:h) {
    temp = (t(x[H:n,]) %*% x[(H - i):(n - i),]) / n
    L = L + (temp %*% t(temp))
  }
  
  # eigenanalysis
  model = eigs_sym(L, k = R)
  E_vec = model$vectors
  E_val = model$values
  if (is.null(r)) {
    if (is.null(R)) {
      R = min(floor(p / 2), floor(n / 2))
    }
    ratios = E_val[2:(R + 1)] / E_val[1:R]
    r = which.min(ratios)
  }
  V = E_vec[,1:r, drop = FALSE]
  
  # Extract factors and residuals
  f_hat = x %*% V
  fitted.val = f_hat %*% t(V)
  e_hat = x - fitted.val
  
  res = list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat,
             "fitted.val" = fitted.val, "r_hat" = r, "means" = means)
  class(res) = "factor_model"
  
  return (res)
}

#' Predict method for factor_model objects (Internal)
#' 
#' @param model An object of class `factor model`
#' @param newdata 
#' 
#' @keywords internal
#' @export
predict.factor_model = function (model, newdata) {
  V = model$V
  mu = model$means
  if (is.matrix(newdata)) {
    z_temp = t(t(newdata) - mu)
    Factor = z_temp %*% V
    z_hat = Factor %*% t(V)
    x_hat = t(t(z_hat) + mu)
  } else if (is.vector(newdata)) {
    z_temp = newdata - mu
    Factor = c(t(z_temp) %*% V)
    z_hat = V %*% Factor
    x_hat = c(mu + z_hat)
  } else {
    stop("predict.factor_model: newdata must be matrix or vector")
  }
  
  return (x_hat)
}

#' Riemannian factor model
#' 
#' @param mfd a manifold object
#' @param x data
#' @param r number of factors (default: NULL, in which case the ratio estimator
#'          is used to select it)
#' @param h number of lags to use in construction of the autocovariance-based 
#'          auxiliary matrix (default: 1)
#' @param mu empirical Frèchet mean (default: NULL, automatically estimated from
#'           data)
#' @param basis an orthonormal basis at mu (default: NULL, automatically chosen)
#' @param demean whether to demean the data first (default: TRUE)
#' @param R maximum number of factors to consider in the ratio estimator (default:
#'          NULL, in which case p/2 is used, where p is the dimension of data)
#' @param frechet_mean_args parameters to pass to Fréchet mean estimation in 
#'                          `frechet_mean`
#' @param log_args optional parameters to pass to `Log_mfd`
#' 
#' @return A `riem_factor` object containing
#' * `V` estimated loading matrix
#' * `mu` empirical Fréchet mean
#' * `E` an orthonormal basis at the Fréchet mean
#' * `f_hat` estimated factor process
#' * `e_hat` residuals
#' * `r_hat` number of factors used
#' * `means` means of each coordinate (in the tangent space)
#' 
#' @references 
#'   Huang, S.-C., Chen, R., & Chen, Y. (2026). 
#'   A Riemannian Factor Model for Manifold-Valued Time Series. 
#'   arXiv preprint arXiv:2607.28385.
#' @examples 
#' 
#' TBA
#' 
#' @export
rfm = function (mfd, x, r = NULL, h = 1, mu = NULL, basis = NULL, demean = TRUE, 
                R = NULL, frechet_mean_args = list(), log_args = list()) {
  if (is.null(mu)) {
    mu = do.call(frechet_mean, c(list("mfd" = mfd, "x" = x), frechet_mean_args))
  }
  
  # construct orthonormal basis
  if (is.null(basis)) {
    basis = basis(mfd, mu)
  }
  
  # construct log-mapped data
  z = do.call(Log_mfd, c(list("mfd" = mfd, "p" = mu, "q" = x), log_args))
  
  # express z in coordinate
  z_coord = tangent_to_vec(mfd, v = z, p = mu, E = basis)
  
  # factor model
  model = LYB_fm(z_coord, r = r, h = h, demean = demean, R = R)
  
  res = list("V" = model$V, "mu" = mu, "E" = basis, "f_hat" = model$f_hat, 
             "e_hat" = model$e_hat, "r_hat" = model$r_hat, "means" = model$means,
             "mfd" = mfd)
  class(res) = c("riem_factor", "factor_model")
  
  return (res)
}

#' Pseudo-prediction method for factor_model objects
#' 
#' @param model An object of class `factor model`
#' @param newdata 
#' 
#' @keywords internal
#' @export
predict.riem_factor = function (model, newdata, log_args = list(), 
                                exp_args = list(), ...) {
  mu = model$mu
  E = model$E
  means = model$means
  mfd = model$mfd
  
  log_x = do.call(Log_mfd, c(list("mfd" = mfd, "p" = mu, "q" = newdata),
                             log_args))
  log_x_vec = tangent_to_vec(mfd, v = log_x, p = mu, E = E)
  
  z_hat_vec = predict.factor_model(model, log_x_vec)
  z_hat = vec_to_tangent(mfd, z_hat_vec, mu, E)
  res = do.call(Exp_mfd, c(list("mfd" = mfd, "p" = mu, "v" = z_hat),
                           exp_args))
  
  return (res)
}




