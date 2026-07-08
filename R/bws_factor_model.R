################################################################################
#### Riemannian factor model on the Bures-Wasserstein manifold #################
################################################################################
#
# These are BWS-specific modeling/evaluation functions, not generic manifold
# primitives, so they're relocated from BWS_util.R largely unchanged rather
# than folded into the geod/Exp_map/Log_map/frechet_mean generic interface.

#' @export
symmetric_to_vector = function(X) {
  return(X[upper.tri(X, diag = TRUE)])
}

#' @export
vector_to_symmetric = function(v, n) {
  res = matrix(NA, n, n)
  res[upper.tri(res, diag = TRUE)] = v
  res[lower.tri(res)] = t(res)[lower.tri(res)]

  return(res)
}

#' @export
is.spd = function(x, tol = 1e-8) {
  if (!isSymmetric.matrix(x)) {
    return(FALSE)
  }
  lambda = eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (all(lambda > tol)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @export
project_to_SPD = function(A, epsilon = 0) {
  # epsilon is to ensure positive-definiteness (occasionally used)

  eig = eigen(A)
  eig$values[eig$values < epsilon] = epsilon
  A_SPD = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  return(A_SPD)
}

#' @export
LYB_fm = function(x, r, h, demean = TRUE) {
  # x: (n by d) observation matrix
  # r: number of factors
  # h: number of lags to use
  # Estimate the factor model of Lam, Yao, and Bathia

  if (demean) {
    mean = colMeans(x)
    x = t(t(x) - colMeans(x))
  }

  # Compute the auxiliary positive definite matrix
  pd = 0
  H = h + 1
  n = nrow(x)
  for (i in 1:h) {
    temp = t(x[H:n, ]) %*% x[(H - i):(n - i), ] / n
    pd = pd + temp %*% t(temp)
  }

  # Eigenanalysis
  model = eigen(pd)
  Evec = model$vectors
  V = Evec[, 1:r] # (d by r)
  if (r == 1) {
    V = as.matrix(V)
  }
  evals = model$values

  # Extract factors and residuals
  f_hat = t(t(V) %*% t(x)) # (n by r)
  e_hat = x - f_hat %*% t(V)

  # Estimate the number of factors
  ratios = evals[2:r] / evals[1:(r - 1)]
  r_hat = which.min(ratios)

  return(list(
    "V" = V, "f_hat" = f_hat, "e_hat" = e_hat,
    "fitted.val" = f_hat %*% t(V), "r_hat" = r_hat,
    "mean" = mean
  ))
}

#' @export
predict_fm = function(V, mu, new_x) {
  if (is.matrix(new_x)) {
    z_temp = t(t(new_x) - mu)
    Factor = z_temp %*% V
    z_hat = Factor %*% t(V)
    x_hat = t(t(z_hat) + mu)
  } else if (is.vector(new_x)) {
    z_temp = new_x - mu
    Factor = c(t(z_temp) %*% V)
    z_hat = V %*% Factor
    x_hat = mu + z_hat
  } else {
    stop("predict_fm: new_x must be matrix or vector")
  }

  return(x_hat)
}

#' Estimate the RFM in Bures-Wasserstein manifold
#'
#' This function estimates the Riemannian factor model (RFM) with data being
#' symmetric positive definite matrices equipped with the Bures-Wasserstein
#' metric.
#'
#' @param x an $n \times p \times p$ array of data where sample is a $p \times p$ SPD matrix
#' @param r number of factors to extract
#' @param h number of lags used in estimating the factor model
#' @return
#' \describe{
#'  \item{A}{Estimated loading matrix}
#'  \item{f_hat}{Estimated factor process}
#'  \item{E}{A set of orthonormal basis for the tangent space}
#'  \item{E_lyapunov}{A set of orthonormal basis passed to the Lyapunov operator}
#'  \item{mu_hat}{Estimated Frechet mean}
#' }
#' @examples
#' set.seed(1)
#' n = 30; p = 3
#' X = array(NA, dim = c(n, p, p))
#' for (i in 1:n) {
#'   A = matrix(rnorm(p * p), p, p)
#'   X[i, , ] = A %*% t(A) + diag(p)
#' }
#'
#' model = rfm_bws(X, r = 1, h = 3)
#' model$mu_hat          # estimated Frechet mean
#' dim(model$f_hat)       # estimated factor process, n by r
#'
#' # evaluate fraction of variance explained on (here, in-sample) data
#' Frac_Var_bws(X[1:5, , ], model, Euclidean_mean = model$mu_hat)
#' @export
rfm_bws = function (x, r, h = 6, batch_size = NULL, max.iter = 100,
                    mu_hat = NULL) {
  n = dim(x)[1]
  p = dim(x)[2]

  # Estimate mu
  if (is.null(mu_hat)) {
    mu_hat = frechet_mean(manifold_bws(), x, batch_size = batch_size, max.iter = max.iter,
                          tau = 0.5, tol = -1, verbose = FALSE)
  }

  # Construct a set of orthonormal basis
  coord = tan_basis_bws(mu_hat)
  E = coord$E
  E_lyapunov = coord$E_lyapunov

  # construct log-mapped data
  log_x_vec = log_vec_construct(x, mu_hat, E_lyapunov)

  # Estimate the factor model
  model = LYB_fm(log_x_vec, r = r, h = h)

  return(list("A" = model$V, "f_hat" = model$f_hat, "E" = E, "E_lyapunov" = E_lyapunov,
              "mu_hat" = mu_hat, "factor_model" = model,
              "r_hat" = model$r_hat))
}

#' Computes evaluation metric for RFM on BWS
#'
#' @param x_test raw test data
#' @param RFM_model output from rfm_bws
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' @return a length-r vector, the fraction of variance unexplained (or squared prediction error) using 1..r factors
#' @examples
#' set.seed(1)
#' n = 30; p = 3
#' X = array(NA, dim = c(n, p, p))
#' for (i in 1:n) {
#'   A = matrix(rnorm(p * p), p, p)
#'   X[i, , ] = A %*% t(A) + diag(p)
#' }
#' model = rfm_bws(X, r = 1, h = 3)
#' Frac_Var_bws(X[1:5, , ], model, Euclidean_mean = model$mu_hat)
#' @export
Frac_Var_bws = function (x_test, RFM_model, evaluation_type = "BWS", fraction = TRUE,
                         Euclidean_mean, return_predictions = FALSE) {

  mu_hat = RFM_model$mu_hat
  E = RFM_model$E
  E_lyapunov = RFM_model$E_lyapunov
  factor_model = RFM_model$factor_model

  V = factor_model$V
  z_mean = factor_model$mean
  r = dim(V)[2]

  if (length(dim(x_test)) == 3) {
    n = dim(x_test)[1]
    p = dim(x_test)[2]
    x.is.array = TRUE
  } else if (length(dim(x_test)) == 2) {
    n = 1
    p = dim(x_test)[1]
    x.is.array = FALSE
  }

  # construct log-mapped data
  log_x_vec = log_vec_construct(x_test, mu_hat, E_lyapunov)

  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], z_mean, log_x_vec)

    # predict
    if (x.is.array) {
      x_hat = array(NA, dim = c(n, p, p))
      for (m in 1:n) {
        temp = log_to_tangent(z_hat[m,], E)
        x_hat[m,,] = Exp_BWS(temp, mu_hat)
      }
    } else {
      temp = log_to_tangent(z_hat, E)
      x_hat = Exp_BWS(temp, mu_hat)
    }

    # evaluate
    if (evaluation_type == "BWS") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + geod_BWS_core(x_hat[m,,], x_test[m,,])^2
          if (i == 1) {
            total_var = total_var + geod_BWS_core(mu_hat, x_test[m,,])^2
          }
        }
      } else {
        res[i] = geod_BWS_core(x_hat, x_test)^2
      }
    } else if (evaluation_type == "Euclidean") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + norm(x_hat[m,,] - x_test[m,,], type = "F")^2
          if (i == 1) {
            total_var = total_var + norm(Euclidean_mean - x_test[m,,], type = "F")^2
          }
        }
      } else {
        res = norm(x_hat - x_test, type = "F")^2
      }
    } else {
      stop("Frac_Var_bws: unsupported evaluation type")
    }
  }

  if (fraction) {
    res = res / total_var
  }

  if (return_predictions) {
    return (list("res" = res, "xhat" = x_hat))
  }

  return (res)
}

#' Computes evaluation metric for linear factor model
#'
#' @param x_test raw test data
#' @param factor_model output from LYB_fm
#' @param mu_hat BWS Frechet mean, externally given
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' @param epsilon used in projecting predictions to SPD (only when evaluation type == "BWS")
#' @export
Frac_Var_LYB = function (x_test, factor_model, mu_hat, Euclidean_mean,
                         evaluation_type = "BWS",
                         fraction = TRUE, return_predictions = FALSE,
                         epsilon = 1e-6) {
  V = factor_model$V
  z_mean = factor_model$mean
  r = dim(V)[2]

  if (length(dim(x_test)) == 3) {
    n = dim(x_test)[1]
    p = dim(x_test)[2]
    x.is.array = TRUE
  } else if (length(dim(x_test)) == 2) {
    n = 1
    p = dim(x_test)[1]
    x.is.array = FALSE
  }

  # encode SPD as vectors
  if (x.is.array) {
    z = array(NA, dim = c(n, p * (p + 1) / 2))
    for (m in 1:n) {
      z[m,] = symmetric_to_vector(x_test[m,,])
    }
  } else {
    z = symmetric_to_vector(x_test)
  }

  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(as.matrix(V[,1:i]), z_mean, z)

    # predict
    if (x.is.array) {
      x_hat = array(NA, dim = c(n, p, p))
      for (m in 1:n) {
        x_hat[m,,] = vector_to_symmetric(z_hat[m,], p)
      }
    } else {
      x_hat = vector_to_symmetric(z_hat, p)
    }

    # evaluate
    if (evaluation_type == "BWS") {
      if (x.is.array) {
        for (m in 1:n) {
          temp = project_to_SPD(x_hat[m,,], epsilon)
          res[i] = res[i] + (Re(geod_BWS_core(temp, x_test[m,,])))^2
          if (i == 1) {
            total_var = total_var + (Re(geod_BWS_core(mu_hat, x_test[m,,])))^2
          }
        }
      } else {
        res[i] = (Re(geod_BWS_core(project_to_SPD(x_hat, epsilon), x_test)))^2
      }
    } else if (evaluation_type == "Euclidean") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + norm(x_hat[m,,] - x_test[m,,], type = "F")^2
          if (i == 1) {
            total_var = total_var + norm(Euclidean_mean - x_test[m,,], type = "F")^2
          }
        }
      } else {
        res = norm(x_hat - x_test, type = "F")^2
      }
    } else {
      stop("Frac_Var_bws: unsupported evaluation type")
    }
  }

  if (fraction) {
    res = res / total_var
  }

  if (return_predictions) {
    return (list("res" = res, "xhat" = x_hat))
  }

  return (res)
}

#' Computes evaluation metric for Oracle on BWS
#'
#' @param x_test raw test data
#' @param dta object output from dta_gen_BWS
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' @export
Frac_Var_ora = function (x_test, dta,
                         evaluation_type = "BWS", fraction = TRUE,
                         Euclidean_mean, return_predictions = FALSE) {

  mu_hat = dta$mu
  E = dta$coord$E
  E_lyapunov = dta$coord$E_lyapunov

  Factors = dta$Factors
  V = dta$A
  r = dim(V)[2]

  if (length(dim(x_test)) == 3) {
    n = dim(x_test)[1]
    p = dim(x_test)[2]
    x.is.array = TRUE
  } else if (length(dim(x_test)) == 2) {
    n = 1
    p = dim(x_test)[1]
    x.is.array = FALSE
  }

  # construct log-mapped data
  log_x_vec = log_vec_construct(x_test, mu_hat, E_lyapunov)

  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], 0, log_x_vec)

    # predict
    if (x.is.array) {
      x_hat = array(NA, dim = c(n, p, p))
      for (m in 1:n) {
        temp = log_to_tangent(z_hat[m,], E)
        x_hat[m,,] = Exp_BWS(temp, mu_hat)
      }
    } else {
      temp = log_to_tangent(z_hat, E)
      x_hat = Exp_BWS(temp, mu_hat)
    }

    # evaluate
    if (evaluation_type == "BWS") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + geod_BWS_core(x_hat[m,,], x_test[m,,])^2
          if (i == 1) {
            total_var = total_var + geod_BWS_core(mu_hat, x_test[m,,])^2
          }
        }
      } else {
        res[i] = geod_BWS_core(x_hat, x_test)^2
      }
    } else if (evaluation_type == "Euclidean") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + norm(x_hat[m,,] - x_test[m,,], type = "F")^2
          if (i == 1) {
            total_var = total_var + norm(Euclidean_mean - x_test[m,,], type = "F")^2
          }
        }
      } else {
        res = norm(x_hat - x_test, type = "F")^2
      }
    } else {
      stop("Frac_Var_bws: unsupported evaluation type")
    }
  }

  if (fraction) {
    res = res / total_var
  }

  if (return_predictions) {
    return (list("res" = res, "xhat" = x_hat))
  }

  return (res)
}

#' @export
VAR1 <- function(X) {
  X <- as.matrix(X)
  k <- ncol(X)
  df <- embed(X, 2)
  Y <- df[, 1:k]
  Z <- cbind(1, df[, -(1:k)])
  B <- solve(t(Z) %*% Z, t(Z) %*% Y)
  return(B)
}

#' Computes predictions of RFM on BWS, with factors predicted by VAR(1)
#'
#' @export
dyn_RFM = function (x, r, test_size = 1, h = 6, batch_size = NULL, max.iter = 100) {

  n = dim(x)[1]
  p = dim(x)[2]
  x_hat = array(NA, dim = c(test_size, p, p))

  for (m in 0:(test_size - 1)) {
    x_train = x[1:(n - test_size + m),,]
    x_test = x[c(1:(n - test_size + m)),,]

    # Get Factors
    if (m == 0) {
      aux = main_BWS(x, r = r, test_size = test_size - m, h = h,
                     batch_size = batch_size, max.iter = max.iter)
      Factors = aux$Factors
      mu_hat = aux$mu_hat
      E = aux$E
      V = aux$V
      z_bar = aux$z_bar
    } else {
      aux = main_BWS(x, r = r, test_size = test_size - m, h = h,
                     batch_size = batch_size, max.iter = max.iter,
                     mu_hat = mu_hat)
      Factors = aux$Factors
      E = aux$E
      V = aux$V
      z_bar = aux$z_bar
    }

    # Predict factors
    B = VAR1(Factors)
    f_hat = as.vector(c(1, tail(Factors, 1)) %*% B)

    # predict
    z_hat = V %*% f_hat + z_bar
    temp = log_to_tangent(z_hat, E)
    x_hat[m + 1,,] = Exp_BWS(temp, mu_hat)

    cat("dyn_RFM: iteration", m, "\n")
  }

  return (x_hat)
}

#' @export
dyn_LFM = function (x, r, test_size = 1, h = 6) {
  n = dim(x)[1]
  p = dim(x)[2]
  x_hat = array(NA, dim = c(test_size, p, p))

  x_vector = array(NA, dim = c(n, p * (p + 1) / 2))
  for (m in 1:n) {
    x_vector[m,] = symmetric_to_vector(x[m,,])
  }

  for (m in 0:(test_size - 1)) {
    x_train = x_vector[1:(n - test_size + m),]

    # Get factors
    aux = LYB_fm(x_train, r = r, h = h)
    Factors = aux$f_hat
    V = aux$V
    zbar = aux$mean


    # predict factors
    B = VAR1(Factors)
    f_hat = as.vector(c(1, tail(Factors, 1)) %*% B)

    # predict
    z_hat = V %*% f_hat + zbar
    x_hat[m + 1,,] = vector_to_symmetric(z_hat, p)
  }

  return (x_hat)
}
