# Essential functions for Bures-Wasserstein factor analysis
library(maotai)
library(expm)
library(deSolve)

geod_BWS_core = function (X, Y) {
    temp = sqrtm(X)
    temp = sqrtm(temp %*% Y %*% temp)
    res = sum(diag(X)) + sum(diag(Y)) - 2 * sum(diag(temp))
    
    return (sqrt(res))
}

geod_BWS = function (X, Y) {
  
  dimX = length(dim(X))
  dimY = length(dim(Y))
  
  if (dimX == 3 && dimY == 3) {
    if (dim(X)[1] != dim(Y)[1]) {
      stop("geod_BWS: X and Y have different sample sizes")
    }
    n = dim(X)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X[i,,], Y[i,,])
    }
  } else if (dimX == 3 && dimY == 2) {
    n = dim(X)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X[i,,], Y)
    }
  } else if (dimX == 2 && dimY == 3) {
    n = dim(Y)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X, Y[i,,])
    }
  } else if (dimX == 2 && dimY == 2) {
    res = geod_BWS_core(X, Y)
  } else {
    stop(paste("geod_BWS: incorrect dimensions X:",
               dim(X),
               "Y:", dim(Y)))
  }
  
  return (res)
}

Exp_BWS_core = function (X, M) {
  L = lyapunov(M ,X)
  p = dim(X)[1]
  diag(L) = diag(L) + 1
  output = L %*% M %*% L
  
  return (output)
}

Exp_BWS = function (X, M) {
  # X: (n by p by p) or (p by p) array of data on the tangent space
  # M: (p by p) reference point
  # tangent space --> sphere
  
  output = array(NA, dim = dim(X))
  
  if (length(dim(X)) == 3) {
    for (i in 1:dim(X)[1]) {
      output[i,,] = Exp_BWS_core(X[i,,], M)
    }
  } else if (length(dim(X)) == 2) {
    output = Exp_BWS_core(X, M)
  } else {
    stop(paste("Exp_BWS: Incorrect dimensions X", dim(X)))
  }
  
  return (output)
}

Log_BWS_core = function (X, M) {
  
  # check if X in the injectivity radius
  test_X = lyapunov(M, X)
  diag(test_X) = diag(test_X) + 1
  eig = eigen(test_X, symmetric = T, only.values = T)$values
  if (!all(eig > 0)) {
    stop("Log_BWS: X outside injectivity radius")
  }
  
  output = sqrtm(X %*% M)
  output = output + sqrtm(M %*% X) - 2 * M
  
  return (output)
}

Log_BWS = function (X, M) {
  if (length(dim(X)) == 3) {
    output = array(NA, dim = dim(X))
    for (i in 1:dim(X)[1]) {
      output[i,,] = Log_BWS_core(X[i,,], M)
    }
  } else if (length(dim(X)) == 2) {
    output = Log_BWS_core(X, M)
  } else {
    stop(paste("Log_BWS: Incorrect dimensions X", dim(X)))
  }
  
  return (output)
}

#' express the log-mapped data in E coordinates
#' 
#' express the tangent vector log_M(x) in a canonical orthonormal basis
log_vec_construct = function (x, M, E_lyapunov = NULL) {
  log_x = Log_BWS(x, M)
  
  if (is.null(E_lyapunov)) {
    E_lyapunov = tan_basis_bws(M)$E_lyapunov
  }
  
  if (length(dim(x)) == 3) {
    n = dim(x)[1]
    p = dim(x)[2]
    log_x_vec = array(NA, dim = c(n, p * (p + 1) / 2))
    
    for (m in 1:n) {
      counter = 0
      for (i in 1:p) {
        for (j in i:p) {
          counter = counter + 1
          log_x_vec[m, counter] = 0.5 * sum(diag(E_lyapunov[counter,,] %*% log_x[m,,]))
        }
      }
    }
  } else if (length(dim(x)) == 2) {
    p = dim(x)[1]
    log_x_vec = rep(NA, p * (p + 1) / 2)
    
    counter = 0
    for (i in 1:p) {
      for (j in i:p) {
        counter = counter + 1
        log_x_vec[counter] = 0.5 * sum(diag(E_lyapunov[counter,,] %*% log_x))
      }
    }
    log_x_vec = t(as.matrix(log_x_vec))
  } else {
    stop(paste("log_vec_construct: Incorrect dimensions", dim(x)))
  }
  
  return (log_x_vec)
}

#' Identify a tangent vector (in E coordinate) as a symmetric matrix
#' 
#' @param z a coordinate in E
#' @return a symmetric matrix in the tangent space
#' 
log_to_tangent = function (z, E) {
  p = dim(E[1,,])[1]
  res = array(0, dim = c(p, p))
  if (length(z) != p * (p + 1) / 2) {
    stop("log_to_tangent: z has wrong length")
  }
  
  counter = 0
  for (i in 1:p) {
    for (j in i:p) {
      counter = counter + 1
      res = res + z[counter] * E[counter,,]
    }
  }
  
  return (res)
}

#' express the tangent vector V (identified as symmetric matrix) in E coordinate
tangent_in_E = function (V, M, E_lyapunov = NULL) {
  p = dim(V)[1]
  if (is.null(E_lyapunov)) {
    E_lyapunov = tan_basis_bws(M)$E_lyapunov
  }
  
  res = rep(NA, p * (p + 1) / 2)
  counter = 0
  for (i in 1:p) {
    for (j in i:p) {
      counter = counter + 1
      res[counter] = 0.5 * sum(diag(E_lyapunov[counter,,] %*% V))
    }
  }
  
  return (res)
}

symmetric_to_vector = function (X) {
  return (X[upper.tri(X, diag = TRUE)])
}

vector_to_symmetric = function (v, n) {
  res = matrix(NA, n, n)
  res[upper.tri(res, diag = TRUE)] = v
  res[lower.tri(res)] = t(res)[lower.tri(res)]
  
  return (res)
}

is.spd = function (x, tol = 1e-8) {
  if (!isSymmetric.matrix(x)) {
    return (FALSE)
  } 
  lambda = eigen(x, symmetric = T, only.values = T)$values
  if (all(lambda > tol)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}

project_to_SPD <- function(A, epsilon = 0) {
  # epsilon is to ensure positive-definiteness (occassionally used)
  
  eig <- eigen(A)
  eig$values[eig$values < epsilon] <- epsilon
  A_SPD <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  return(A_SPD)
}

mean_on_BWS = function (X, tau = 0.1, tol = 1e-4, max.iter = 1000,
                        batch_size = NULL, verbose = FALSE) {
  # X: (n by ...) data
  # tau: step size
  tau_0 = tau
  
  n = dim(X)[1]
  p = dim(X)[2]
  mu = X[1,,]
  for (i in 1:max.iter) {
    # mini-batch sampling
    if (!is.null(batch_size)) {
      idx = sample(n, batch_size, replace = FALSE)
      X_batch = X[idx,,]
      grad = apply(Log_BWS(X_batch, mu), c(2, 3), mean) # stochastic gradient
    } else {
      grad = colMeans(matrix(Log_BWS(X, mu), nrow = dim(X)[1]))
      grad = matrix(grad, nrow = p, ncol = p)
    }
    # grad = apply(Log_BWS(X, mu), c(2, 3), mean)
    mu_new = Exp_BWS(tau * grad, mu)
    
    loss = mean(geod_BWS(X, mu_new))
    if (verbose) {
      cat("mean_on_BWS: iter", i, "loss", round(loss, 4), "\n")
    }
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
    
    if (!is.null(batch_size)) {
      tau = tau_0 / sqrt(i)
    }
  }
  
  return (mu)
}

LYB_fm = function (x, r, h, demean = TRUE) {
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
    temp = t(x[H:n,]) %*% x[(H - i):(n - i),] / n
    pd = pd + temp %*% t(temp)
  }
  
  # Eigenanalysis
  model = eigen(pd)
  Evec = model$vectors
  V = Evec[,1:r] # (d by r)
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
  
  return (list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat, 
               "fitted.val" = f_hat %*% t(V), "r_hat" = r_hat,
               "mean" = mean))
}

predict_fm = function (V, mu, new_x) {
  
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
  
  return (x_hat)
}

Christoffel_BWS_core = function (Sigma, X, Y) {
  Lx = lyapunov(Sigma, X)
  Ly = lyapunov(Sigma, Y)
  
  res = (Sigma %*% Ly %*% Lx) + (Ly %*% Lx %*% Sigma) - (Lx %*% Y) - (Ly %*% X)
  res = 0.5 * (res + t(res))
  
  return (res)
}

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
#' @export
pt_bws = function (Sigma1, Sigma2, V, method = "adams") {
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

#' Computes a canonical basis for the tangent space at a point in BWS
#' 
#' @param Sigma Base point
#' @return An array of orthonormal basis (E), and an array of E passed through the Lyapunov operator
#' @export
tan_basis_bws = function (Sigma) {
  p = dim(Sigma)[1]
  
  model = eigen(Sigma)
  lambdas = model$values
  P = model$vectors
  
  # construct canonical basis
  E = array(NA, dim = c(p * (p + 1) / 2, p, p))
  E_lyapunov = array(NA, dim = c(p * (p + 1) / 2, p, p))
  counter = 0
  for (i in 1:p) {
    for (j in i:p) {
      counter = counter + 1
      S = matrix(0, ncol = p, nrow = p)
      S_tilde = matrix(0, ncol = p, nrow = p)
      if (i == j) {
        S[i,j] = sqrt(2 * (lambdas[i] + lambdas[j]))
        S_tilde[i,j] = 1 / sqrt(lambdas[i])
      } else {
        S[i,j] = sqrt(lambdas[i] + lambdas[j])
        S[j,i] = sqrt(lambdas[i] + lambdas[j])
        
        S_tilde[i,j] = 1 / sqrt(lambdas[i] + lambdas[j])
        S_tilde[j,i] = 1 / sqrt(lambdas[i] + lambdas[j])
      }
      E[counter,,] = P %*% S %*% t(P)
      E_lyapunov[counter,,] = P %*% S_tilde %*% t(P) # L_{mu_hat}(E)
    }
  }
  
  return (list("E" = E, "E_lyapunov" = E_lyapunov))
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
#' @param mu_tol optimization tolerance level for estimating the Fr\'{e}chet mean
#' @return 
#' \describe{
#'  \item{A}{Estimated loading matrix}
#'  \item{f_hat}{Estimated factor process}
#'  \item{E}{A set of orthonormal basis for the tangent space}
#'  \item{E_lyapunov}{A set of orthonormal basis passed to the Lyapunov operator}
#'  \item{mu_hat}{Estimated Fr\'{e}chet mean}
#' } 
#' @export
rfm_bws = function (x, r, h = 6, batch_size = NULL, max.iter = 100,
                    mu_hat = NULL) {
  n = dim(x)[1]
  p = dim(x)[2]
  
  # Estimate mu
  if (is.null(mu_hat)) {
    mu_hat = mean_on_BWS(x, batch_size = batch_size, max.iter = max.iter,
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


#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for RFM on BWS
#' 
#' @param x_test raw test data
#' @param factor_model output from rfm_bws; an LYB_fm object
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' 
Frac_Var_bws = function (x_test, RFM_model, evaluation_type = "BWS", fraction = TRUE,
                         Euclidean_mean, return_predictions = FALSE) {
  
  mu_hat = RFM_model$mu_hat
  E = RFM_model$E
  E_lyapunov = RFM_model$E_lyapunov
  factor_model = RFM_model$factor_model
  Factors = RFM_model$f_hat
  
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

#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for linear factor model
#' 
#' @param x_test raw test data
#' @param factor_model output from LYB_fm
#' @param mu_hat BWS Frechet mean, externally given
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' @param epsilon used in projecting predictions to SPD (only when evaluation type == "BWS")
#' 
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

#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for Oracle on BWS
#' 
#' @param x_test raw test data
#' @param dta object output from dta_gen_BWS
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' 
Frac_Var_ora = function (x_test, dta,
                         evaluation_type = "BWS", fraction = TRUE,
                         Euclidean_mean, return_predictions = FALSE) {
  
  mu_hat = dta$mu
  E = dta$coord$E
  E_lyapunov = dta$coord$E_lyapunov
  
  # factor_model = RFM_model$factor_model
  
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


VAR1 <- function(X) {
  X <- as.matrix(X)
  k <- ncol(X)
  df <- embed(X, 2)                    
  Y <- df[, 1:k]                       
  Z <- cbind(1, df[, -(1:k)])          
  B <- solve(t(Z) %*% Z, t(Z) %*% Y)   
  return(B)
}

#' This function is used internally to evaluate the factor models
#' 
#' Computes predictions of RFM on BWS, with factors predicted by VAR(1)
#' 
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



