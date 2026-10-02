# Generic Fréchet mean via Riemannian Gradient Descent

determine_batch = function (n, batch_size) {
  if (batch_size <= 1) {
    return (max(1, round(batch_size * n)))
  }
  return (min(n, round(batch_size)))
}

is_batch = function (mfd, v) {
  if (mfd$point_ndim == 1) {
    if (!is.matrix(x)) {
      return (FALSE)
    } else if (nrow(x) == 1) {
      return (FALSE)
    } else {
      return (TRUE)
    }
  } else {
    if (length(dim(v)) == mfd$point_ndim) {
      return (FALSE)
    } else if (dim(v)[1] == 1) {
      return (FALSE)
    } else {
      return (TRUE)
    }
  }
}

# average a batch of tangent vectors/matrices produced by Log_map()
mean_tangent = function (mfd, v) {
  if (!is_batch(mfd, v)) {
    return(v)
  }
  if (mfd$point_ndim == 1) {
    return(colMeans(v))
  } else {
    return(apply(v, seq_len(mfd$point_ndim) + 1, mean))
  }
}


#' Fréchet mean on a Riemannian manifold
#' 
#' A single generic implementation shared by all manifolds: at each iteration,
#' average the Log-mapped tangent vectors of a (sub)sample of the data at the
#' current estimate, take an Exp step in that direction, and stop once the drop
#' in average geodesic distance falls below tolerance
#' 
#' @param mfd a manifold object
#' @param x an (n by ...) data
#' @param tau step size
#' @param tol convergence tolerance
#' @param max.iter maximum number of iterations
#' @param batch_size NULL for full-batch gradient descent; a value in (0, 1]
#'    for a proportion of the data; a value > 1 for an absolute subsample count
#'    (stochastic gradient descent, with tau decaying as tau / sqrt(iter))
#' @param init optional initial value; if NULL a random data point is used
#' @param verbose if TRUE, print the loss at each iteration
#' 
#' @return the estimated Fréchet mean, a single point on the manifold
#' 
#' @export
frechet_mean.default = function (mfd, x, tau = 0.1, tol = 1e-6, max.iter = 1000,
                                 batch_size = NULL, init = NULL, verbose = FALSE) {
  n = n_points(mfd, x)
  if (n == 1) {
    return (get_point(mfd, x, 1))
  }
  
  if (is.null(init)) {
    mu = get_point(mfd, x, sample(n, 1))
  } else {
    mu = init
  }
  
  tau_0 = tau
  use_batch = !is.null(batch_size)
  
  for (i in 1:max.iter) {
    if (use_batch) {
      bn = determine_batch(n, batch_size)
      idx = sample(n, bn, replace = FALSE)
      x_batch = subset_points(mfd, x, idx)
    } else {
      x_batch = x
    }
    
    grad = mean_tangent(mfd, Log_mfd(mfd, mu, x_batch))
    mu_new = Exp_mfd(mfd, mu, tau * grad)
    
    loss = mean(geod(mfd, x, mu_new)^2)
    if (verbose) {
      cat("frechet_mean:", class(mfd)[1], "iter", i, "loss",
          sprintf("%.4f", loss), "\n")
    }
    
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
    
    if (use_batch) {
      tau = tau_0 / sqrt(i)
    }
  }
  
  return (mu)
}

#' Fréchet mean on a Riemannian manifold
#' 
#' A generic 
#' 
#' @param mfd a manifold object
#' @param x an (n by ...) data
#' @param method either "specialized" (default), which uses specialized algorithm
#'               (if supported), or "SGD" which uses generic Riemannian SGD
#' @param tau step size
#' @param tol covergence tolerance
#' @param max.iter maximum number of iterations
#' @param batch_size NULL for full-batch gradient descent; a value in (0, 1]
#'    for a proportion of the data; a value > 1 for an absolute subsample count
#'    (stochastic gradient descent, with tau decaying as tau / sqrt(iter))
#' @param init optional initial value; if NULL a random data point is used
#' @param verbose if TRUE, print the loss at each iteration
#' 
#' @return the estimated Fréchet mean, a single point on the manifold
#' 
#' @examples 
#' mfd = manifold_bws()
#' x = array(NA, dim = c(10, 3, 3))
#' for (i in 1:10) x[i,,] = crossprod(matrix(rnorm(9), ncol = 3)) + diag(0.1, 3)
#' frechet_mean(mfd, x, verbose = TRUE)
#' frechet_mean(mfd, x, method = "SGD", verbose = TRUE)
#' 
#' @export
frechet_mean = function (mfd, x, ...) {
  UseMethod("frechet_mean")
}






















