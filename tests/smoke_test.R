#!/usr/bin/env Rscript
# Dependency-free smoke test for the manifoldstats generic interface.
# No testthat: plain stopifnot() assertions, run via `Rscript tests/smoke_test.R`.
#
# For each manifold, checks:
#   1. Exp_map(Log_map(x, mu), mu) round-trips back to x
#   2. geod(x, x) == 0
#   3. frechet_mean() matches a hand-copied reference of the *original*
#      pre-refactor mean_on_BWS/mean_on_gr/mean_on_sphere loop structure
#      (same init, no batching, so there's no RNG-driven divergence)
#   4. parallel_transport() runs and returns a same-shaped tangent vector

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
X_bws = array(NA, dim = c(n, p, p))
for (i in 1:n) X_bws[i, , ] = rand_spd(p)
mfd_bws = manifold_bws()

v = Log_map(mfd_bws, X_bws[1, , ], X_bws[2, , ]) # tangent at X_bws[2,,], for the round trip below
x_rt = Exp_map(mfd_bws, v, X_bws[2, , ])
check("BWS round-trip Exp(Log(x, mu), mu) == x", max(abs(x_rt - X_bws[1, , ])) < 1e-6)

# sqrtm() is iterative, so geod(x,x) is only zero up to solver tolerance
check("BWS geod(x, x) == 0", geod(mfd_bws, X_bws[1, , ], X_bws[1, , ]) < 1e-5)

v_pt = Log_map(mfd_bws, X_bws[2, , ], X_bws[1, , ]) # tangent at X_bws[1,,] ("from"), for parallel transport below

mean_on_BWS_orig = function(X, tau, tol, max.iter, init) {
  p = dim(X)[2]
  mu = init
  for (i in 1:max.iter) {
    grad = colMeans(matrix(Log_BWS(X, mu), nrow = dim(X)[1]))
    grad = matrix(grad, nrow = p, ncol = p)
    mu_new = Exp_BWS(tau * grad, mu)
    loss = mean(geod_BWS(X, mu_new))
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    mu = mu_new
    loss_old = loss
  }
  return(mu)
}

mu_ref = mean_on_BWS_orig(X_bws, tau = 0.3, tol = 1e-7, max.iter = 300, init = X_bws[1, , ])
mu_new = frechet_mean(mfd_bws, X_bws, tau = 0.3, tol = 1e-7, max.iter = 300, init = X_bws[1, , ])
check("BWS frechet_mean matches pre-refactor reference", max(abs(mu_ref - mu_new)) < 1e-6)

pt_res = parallel_transport(mfd_bws, X_bws[1, , ], X_bws[2, , ], v_pt)
check("BWS parallel_transport returns p x p matrix", all(dim(pt_res) == c(p, p)))

# backward-compat wrappers agree with the new generic
check("BWS compat geod_BWS matches geod()", abs(geod_BWS(X_bws[1, , ], X_bws[2, , ]) - geod(mfd_bws, X_bws[1, , ], X_bws[2, , ])) < 1e-10)

################################################################################
## Grassmannian
################################################################################

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

v = Log_map(mfd_sp, X_sp[1, ], X_sp[2, ]) # tangent at X_sp[2,], for the round trip below
x_rt = Exp_map(mfd_sp, v, X_sp[2, ])
check("Sphere round-trip Exp(Log(x, mu), mu) == x", max(abs(x_rt - X_sp[1, ])) < 1e-6)

check("Sphere geod(x, x) == 0", geod(mfd_sp, X_sp[1, ], X_sp[1, ]) < 1e-8)

v_pt = Log_map(mfd_sp, X_sp[2, ], X_sp[1, ]) # tangent at X_sp[1,] ("from"), for parallel transport below

mean_on_sphere_orig = function(x, tau, tol, max.iter, init) {
  mu = init
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(tau * grad, mu)
    loss = mean(geod_sphere(x, mu_new))
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
check("Sphere frechet_mean matches pre-refactor reference", max(abs(mu_ref - mu_new)) < 1e-6)

pt_res = parallel_transport(mfd_sp, X_sp[1, ], X_sp[2, ], v_pt)
check("Sphere parallel_transport returns length-q vector", length(pt_res) == q)

################################################################################

if (fail) {
  cat("\nSMOKE TEST FAILED\n")
  quit(status = 1)
} else {
  cat("\nAll smoke tests passed.\n")
}
