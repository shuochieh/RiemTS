################################################################################
#### Backward-compatible wrappers for the old per-manifold function names #####
################################################################################
#
# Nothing in this project currently calls these old names (verified against
# the project's own backup zip), but they're kept as one-line wrappers around
# the new generic interface as a safety net for any external callers not
# visible here. All the real logic now lives in geod()/Exp_map()/Log_map()/
# frechet_mean()/parallel_transport() plus the per-manifold R/*.R files.

# --- Bures-Wasserstein ------------------------------------------------------

#' @export
geod_BWS = function(X, Y) geod(manifold_bws(), X, Y)

#' @export
Exp_BWS = function(X, M) Exp_map(manifold_bws(), X, M)

#' @export
Log_BWS = function(X, M) Log_map(manifold_bws(), X, M)

#' @export
log_to_tangent = function(z, E) tangent_from_coords(E, z)

#' @export
tangent_in_E = function(V, M, E_lyapunov = NULL) {
  if (is.null(E_lyapunov)) {
    E_lyapunov = tan_basis_bws(M)$E_lyapunov
  }
  coords_from_basis(E_lyapunov, V)
}

#' @export
mean_on_BWS = function(X, tau = 0.1, tol = 1e-4, max.iter = 1000,
                       batch_size = NULL, verbose = FALSE) {
  frechet_mean(manifold_bws(), X, tau = tau, tol = tol, max.iter = max.iter,
              batch_size = batch_size, verbose = verbose)
}

#' @export
pt_bws = function(Sigma1, Sigma2, V, method = "adams") {
  parallel_transport(manifold_bws(), Sigma1, Sigma2, V, method = method)
}

# --- Grassmannian ------------------------------------------------------------

#' @export
geod_gr = function(x, y, p = NULL) geod(manifold_grassmann(p), x, y)

#' @export
Exp_gr = function(D, P) Exp_map(manifold_grassmann(), D, P)

#' @export
Log_gr = function(P, x, p = NULL) Log_map(manifold_grassmann(p), x, P)

#' @export
coord_gr = function(basis, x) coords_from_basis(basis, x)

#' @export
coord2tan = function(basis, coord) tangent_from_coords(basis, coord)

#' @export
mean_on_gr = function(x, p, tau = 0.1, tol = 1e-8, max.iter = 1000,
                      batch_por = 1.0, init = NULL, verbose = FALSE) {
  frechet_mean(manifold_grassmann(p), x, tau = tau, tol = tol, max.iter = max.iter,
              batch_size = batch_por, init = init, verbose = verbose)
}

# --- Sphere --------------------------------------------------------------

#' @export
geod_sphere = function(x, y) geod(manifold_sphere(), x, y)

#' @export
Exp_sphere = function(x, mu) Exp_map(manifold_sphere(), x, mu)

#' @export
Log_sphere = function(x, mu) Log_map(manifold_sphere(), x, mu)

#' @export
mean_on_sphere = function(x, tau = 0.1, tol = 1e-8, max.iter = 1000,
                          init = NULL, verbose = FALSE) {
  frechet_mean(manifold_sphere(), x, tau = tau, tol = tol, max.iter = max.iter,
              init = init, verbose = verbose)
}

#' @export
pt_sphere = function(x, y, V) parallel_transport(manifold_sphere(), x, y, V)
