#!/usr/bin/env Rscript
# Dependency-free smoke test for the RiemTS generic manifold interface.
# No testthat: plain stopifnot() assertions, run via `Rscript tests/smoke_test.R`.
#
# For each manifold, checks:
#   1. Exp_mfd(mu, Log_mfd(mu, x)) round-trips back to x
#   2. geod(x, x) == 0
#   3. the generic gradient-descent frechet_mean() matches a hand-written
#      reference loop (same init, no batching, so there's no RNG-driven divergence)
#   4. ptransport() runs and returns a same-shaped tangent vector

pkgload::load_all(".", quiet = TRUE)
set.seed(20260708)

fail = FALSE
check = function(label, cond) {
  status = if (isTRUE(cond)) "OK" else "FAIL"
  if (!isTRUE(cond)) fail <<- TRUE
  cat(sprintf("[%s] %s\n", status, label))
}

################################################################################
## Bures-Wasserstein
################################################################################

rand_spd = function(p) {
  A = matrix(rnorm(p * p), p, p)
  A %*% t(A) + diag(p) * 0.5
}

p = 4
n = 15
X_bw = array(NA, dim = c(n, p, p))
for (i in 1:n) X_bw[i, , ] = rand_spd(p)
mfd_bw = manifold_BW()

v = Log_mfd(mfd_bw, X_bw[2, , ], X_bw[1, , ]) # tangent at X_bw[2,,], for the round trip below
x_rt = Exp_mfd(mfd_bw, X_bw[2, , ], v)
check("BW round-trip Exp(Log(x, mu), mu) == x", max(abs(x_rt - X_bw[1, , ])) < 1e-6)

# sqrtm() is iterative, so geod(x,x) is only zero up to solver tolerance
check("BW geod(x, x) == 0", geod(mfd_bw, X_bw[1, , ], X_bw[1, , ]) < 1e-5)

v_pt = Log_mfd(mfd_bw, X_bw[1, , ], X_bw[2, , ]) # tangent at X_bw[1,,] ("from"), for parallel transport below

mean_on_BW_orig = function(X, tau, tol, max.iter, init) {
  p = dim(X)[2]
  mu = init
  for (i in 1:max.iter) {
    grad = colMeans(matrix(Log_BW(mu, X), nrow = dim(X)[1]))
    grad = matrix(grad, nrow = p, ncol = p)
    mu_new = Exp_BW(tau * grad, mu)
    loss = mean(geod_BW(X, mu_new)^2)
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
  }
  return(mu)
}

mu_ref = mean_on_BW_orig(X_bw, tau = 0.3, tol = 1e-7, max.iter = 300, init = X_bw[1, , ])
mu_new = frechet_mean(mfd_bw, X_bw, method = "SGD", tau = 0.3, tol = 1e-7, max.iter = 300, init = X_bw[1, , ])
check("BW frechet_mean (SGD) matches reference loop", max(abs(mu_ref - mu_new)) < 1e-6)

# the specialized fixed-point algorithm should land on the same barycenter
mu_fp = frechet_mean(mfd_bw, X_bw, tol = 1e-12, max.iter = 1000)
check("BW frechet_mean (fixed point) agrees with SGD", max(abs(mu_fp - mu_new)) < 1e-3)

pt_res = ptransport(mfd_bw, X_bw[1, , ], X_bw[2, , ], v_pt)
check("BW ptransport returns p x p matrix", all(dim(pt_res) == c(p, p)))

# backward-compat wrappers agree with the new generic
check("BW compat geod_BW matches geod()", abs(geod_BW(X_bw[1, , ], X_bw[2, , ]) - geod(mfd_bw, X_bw[1, , ], X_bw[2, , ])) < 1e-10)

################################################################################
## Grassmannian
################################################################################

# Skipped: the Grassmannian has not been ported to the Exp_mfd/Log_mfd/geod/
# ptransport generics yet (see "Manifold operations for the Grassmannian" in the
# README roadmap). Re-enable these checks once R/grassmann.R implements them.
if (FALSE) {
  rand_projector = function(d, p) {
    A = matrix(rnorm(d * d), d, d)
    Q = qr.Q(qr(A))
    U = Q[, 1:p, drop = FALSE]
    U %*% t(U)
  }

  d = 6
  p_gr = 2
  n_gr = 15
  X_gr = array(NA, dim = c(n_gr, d, d))
  for (i in 1:n_gr) X_gr[i, , ] = rand_projector(d, p_gr)
  mfd_gr = manifold_grassmann(p_gr)

  v = Log_map(mfd_gr, X_gr[1, , ], X_gr[2, , ]) # tangent at X_gr[2,,], for the round trip below
  x_rt = Exp_map(mfd_gr, v, X_gr[2, , ])
  check("Grassmann round-trip Exp(Log(x, mu), mu) == x", max(abs(x_rt - X_gr[1, , ])) < 1e-6)

  check("Grassmann geod(x, x) == 0", geod(mfd_gr, X_gr[1, , ], X_gr[1, , ]) < 1e-5)

  v_pt = Log_map(mfd_gr, X_gr[2, , ], X_gr[1, , ]) # tangent at X_gr[1,,] ("from"), for parallel transport below

  mean_on_gr_orig = function(x, p, tau, tol, max.iter, init) {
    mu = init
    for (i in 1:max.iter) {
      basis = basis_gr(mu, p = p)
      temp = Log_gr(mu, x, p = p)
      coord = coord_gr(basis, temp)
      grad_coord = colMeans(coord)
      grad = coord2tan(basis, grad_coord)
      mu_new = Exp_gr(tau * grad, mu)
      loss = mean(geod_gr(x, mu_new, p = p))
      if (i > 1 && (loss_old - loss < tol)) {
        mu = mu_new
        break
      }
      mu = mu_new
      loss_old = loss
    }
    return(mu)
  }

  mu_ref = mean_on_gr_orig(X_gr, p_gr, tau = 0.3, tol = 1e-9, max.iter = 300, init = X_gr[1, , ])
  mu_new = frechet_mean(mfd_gr, X_gr, tau = 0.3, tol = 1e-9, max.iter = 300, init = X_gr[1, , ])
  check("Grassmann frechet_mean matches pre-refactor reference (validates dropping the coord/basis detour)",
        max(abs(mu_ref - mu_new)) < 1e-6)

  pt_res = parallel_transport(mfd_gr, X_gr[1, , ], X_gr[2, , ], v_pt)
  check("Grassmann parallel_transport returns d x d matrix", all(dim(pt_res) == c(d, d)))
}

################################################################################
## Sphere
################################################################################

rand_sphere = function(q) {
  v = rnorm(q)
  v / sqrt(sum(v^2))
}

q = 5
n_sp = 15
X_sp = matrix(NA, n_sp, q)
for (i in 1:n_sp) X_sp[i, ] = rand_sphere(q)
mfd_sp = manifold_sphere()

v = Log_mfd(mfd_sp, X_sp[2, ], X_sp[1, ]) # tangent at X_sp[2,], for the round trip below
x_rt = Exp_mfd(mfd_sp, X_sp[2, ], v)
check("Sphere round-trip Exp(Log(x, mu), mu) == x", max(abs(x_rt - X_sp[1, ])) < 1e-6)

check("Sphere geod(x, x) == 0", geod(mfd_sp, X_sp[1, ], X_sp[1, ]) < 1e-8)

v_pt = Log_mfd(mfd_sp, X_sp[1, ], X_sp[2, ]) # tangent at X_sp[1,] ("from"), for parallel transport below

mean_on_sphere_orig = function(x, tau, tol, max.iter, init) {
  mu = init
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(tau * grad, mu)
    loss = mean(geod_sphere(x, mu_new)^2)
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
  }
  return(mu)
}

mu_ref = mean_on_sphere_orig(X_sp, tau = 0.3, tol = 1e-9, max.iter = 300, init = X_sp[1, ])
mu_new = frechet_mean(mfd_sp, X_sp, tau = 0.3, tol = 1e-9, max.iter = 300, init = X_sp[1, ])
check("Sphere frechet_mean matches reference loop", max(abs(mu_ref - mu_new)) < 1e-6)

pt_res = ptransport(mfd_sp, X_sp[1, ], X_sp[2, ], v_pt)
check("Sphere ptransport returns length-q vector", length(pt_res) == q)

################################################################################

if (fail) {
  cat("\nSMOKE TEST FAILED\n")
  quit(status = 1)
} else {
  cat("\nAll smoke tests passed.\n")
}
