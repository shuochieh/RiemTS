################################################################################
#### Generic Frechet mean via Riemannian gradient descent ######################
################################################################################

# batch_size <= 1 is interpreted as a proportion of n (matches the old
# Grassmannian batch_por convention, where the default 1.0 meant "full
# batch"); batch_size > 1 is interpreted as an absolute sample count
# (matches the old BWS batch_size convention).
determine_batch_n = function(n, batch_size) {
  if (batch_size <= 1) {
    return(max(1, round(batch_size * n)))
  }
  return(min(n, round(batch_size)))
}

#' Frechet mean on a manifold via Riemannian (stochastic) gradient descent
#'
#' A single generic implementation shared by all manifolds: at each
#' iteration, average the Log-mapped tangent vectors of a (sub)sample of the
#' data at the current estimate, take an Exp step in that direction, and
#' stop once the drop in average geodesic distance falls below `tol`.
#'
#' @param mfd a manifold object
#' @param x an (n by ...) batch of data
#' @param tau step size
#' @param tol convergence tolerance on the loss decrease
#' @param max.iter maximum number of iterations
#' @param batch_size NULL for full-batch gradient descent; a value in (0, 1]
#'   for a proportion of the data; a value > 1 for an absolute subsample
#'   count (stochastic gradient descent, with tau decaying as tau / sqrt(iter))
#' @param init optional initial value; if NULL a random data point is used
#' @param verbose if TRUE, print the loss at each iteration
#' @export
frechet_mean = function(mfd, x, tau = 0.1, tol = 1e-8, max.iter = 1000,
                        batch_size = NULL, init = NULL, verbose = FALSE) {
  n = n_points(mfd, x)
  if (n == 1) {
    return(get_point(mfd, x, 1))
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
      bn = determine_batch_n(n, batch_size)
      idx = sample(n, bn, replace = FALSE)
      x_batch = subset_points(mfd, x, idx)
    } else {
      x_batch = x
    }

    grad = mean_tangent(mfd, Log_map(mfd, x_batch, mu))
    mu_new = Exp_map(mfd, tau * grad, mu)

    loss = mean(geod(mfd, x, mu_new))
    if (verbose) {
      cat("frechet_mean:", class(mfd)[1], "iter", i, "loss", round(loss, 6), "\n")
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

  return(mu)
}
