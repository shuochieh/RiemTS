################################################################################
#### Backward-compatible wrappers for the old per-manifold function names #####
################################################################################
#
# Nothing in this project currently calls these old names (verified against
# the project's own backup zip), but they're kept as one-line wrappers around
# the new generic interface as a safety net for any external callers not
# visible here. All the real logic now lives in geod()/Exp_map()/Log_map()/
# frechet_mean()/parallel_transport() plus the per-manifold R/*.R files. New
# code should call those directly; these are documented here only to point
# at their replacement.

# --- Bures-Wasserstein ------------------------------------------------------

#' Legacy alias for `geod(manifold_bws(), X, Y)`
#' @param X,Y single SPD matrices or (n by p by p) batches
#' @seealso [geod()]
#' @export
geod_BWS = function(X, Y) geod(manifold_bws(), X, Y)

#' Legacy alias for `Exp_map(manifold_bws(), X, M)`
#' @param X tangent vector(s) at M
#' @param M base point
#' @seealso [Exp_map()]
#' @export
Exp_BWS = function(X, M) Exp_map(manifold_bws(), X, M)

#' Legacy alias for `Log_map(manifold_bws(), X, M)`
#' @param X point(s) on the manifold
#' @param M base point
#' @seealso [Log_map()]
#' @export
Log_BWS = function(X, M) Log_map(manifold_bws(), X, M)

#' Legacy alias for `tangent_from_coords(E, z)`
#' @param z coordinate vector
#' @param E basis array
#' @seealso [tangent_from_coords()]
#' @export
log_to_tangent = function(z, E) tangent_from_coords(E, z)

#' Legacy alias for `coords_from_basis(E_lyapunov, V)`
#' @param V a tangent vector (symmetric matrix), or a batch of them
#' @param M base point, used to compute E_lyapunov if not supplied
#' @param E_lyapunov optional precomputed basis from `tan_basis_bws(M)$E_lyapunov`
#' @seealso [coords_from_basis()]
#' @export
tangent_in_E = function(V, M, E_lyapunov = NULL) {
  if (is.null(E_lyapunov)) {
    E_lyapunov = tan_basis_bws(M)$E_lyapunov
  }
  coords_from_basis(E_lyapunov, V)
}

#' Legacy alias for `frechet_mean(manifold_bws(), X, ...)`
#' @param X an (n by p by p) array of data
#' @param tau step size
#' @param tol convergence tolerance
#' @param max.iter maximum number of iterations
#' @param batch_size NULL for full-batch, or an absolute subsample count
#' @param verbose if TRUE, print the loss at each iteration
#' @seealso [frechet_mean()]
#' @export
mean_on_BWS = function(X, tau = 0.1, tol = 1e-4, max.iter = 1000,
                       batch_size = NULL, verbose = FALSE) {
  frechet_mean(manifold_bws(), X, tau = tau, tol = tol, max.iter = max.iter,
              batch_size = batch_size, verbose = verbose)
}

#' Legacy alias for `parallel_transport(manifold_bws(), Sigma1, Sigma2, V, method)`
#' @param Sigma1 starting point
#' @param Sigma2 end point
#' @param V tangent vector at Sigma1
#' @param method ODE solver method passed to `deSolve::ode`
#' @seealso [parallel_transport()]
#' @export
pt_bws = function(Sigma1, Sigma2, V, method = "adams") {
  parallel_transport(manifold_bws(), Sigma1, Sigma2, V, method = method)
}

# --- Grassmannian ------------------------------------------------------------

#' Legacy alias for `geod(manifold_grassmann(p), x, y)`
#' @param x,y single projectors or (n by d by d) batches
#' @param p rank of the projectors; inferred from data if NULL
#' @seealso [geod()]
#' @export
geod_gr = function(x, y, p = NULL) geod(manifold_grassmann(p), x, y)

#' Legacy alias for `Exp_map(manifold_grassmann(), D, P)`
#' @param D tangent vector(s) at P
#' @param P base projector
#' @seealso [Exp_map()]
#' @export
Exp_gr = function(D, P) Exp_map(manifold_grassmann(), D, P)

#' Legacy alias for `Log_map(manifold_grassmann(p), x, P)`
#' @param P base projector
#' @param x projector(s) on the manifold
#' @param p rank of the projectors; inferred from data if NULL
#' @seealso [Log_map()]
#' @export
Log_gr = function(P, x, p = NULL) Log_map(manifold_grassmann(p), x, P)

#' Legacy alias for `coords_from_basis(basis, x)`
#' @param basis an (m by d by d) array of basis matrices
#' @param x a (d by d) matrix or (n by d by d) array
#' @seealso [coords_from_basis()]
#' @export
coord_gr = function(basis, x) coords_from_basis(basis, x)

#' Legacy alias for `tangent_from_coords(basis, coord)`
#' @param basis an (m by d by d) array of basis matrices
#' @param coord coordinate vector, or (n by m) matrix
#' @seealso [tangent_from_coords()]
#' @export
coord2tan = function(basis, coord) tangent_from_coords(basis, coord)

#' Legacy alias for `frechet_mean(manifold_grassmann(p), x, ...)`
#' @param x an (n by d by d) array of data
#' @param p rank of the projectors
#' @param tau step size
#' @param tol convergence tolerance
#' @param max.iter maximum number of iterations
#' @param batch_por proportion of data subsampled per iteration (1.0 = full batch)
#' @param init optional initial value
#' @param verbose if TRUE, print the loss at each iteration
#' @seealso [frechet_mean()]
#' @export
mean_on_gr = function(x, p, tau = 0.1, tol = 1e-8, max.iter = 1000,
                      batch_por = 1.0, init = NULL, verbose = FALSE) {
  frechet_mean(manifold_grassmann(p), x, tau = tau, tol = tol, max.iter = max.iter,
              batch_size = batch_por, init = init, verbose = verbose)
}

# --- Sphere --------------------------------------------------------------

#' Legacy alias for `geod(manifold_sphere(), x, y)`
#' @param x,y single points or (n by q) batches on the sphere
#' @seealso [geod()]
#' @export
geod_sphere = function(x, y) geod(manifold_sphere(), x, y)

#' Legacy alias for `Exp_map(manifold_sphere(), x, mu)`
#' @param x tangent vector(s) at mu
#' @param mu base point
#' @seealso [Exp_map()]
#' @export
Exp_sphere = function(x, mu) Exp_map(manifold_sphere(), x, mu)

#' Legacy alias for `Log_map(manifold_sphere(), x, mu)`
#' @param x point(s) on the sphere
#' @param mu base point
#' @seealso [Log_map()]
#' @export
Log_sphere = function(x, mu) Log_map(manifold_sphere(), x, mu)

#' Legacy alias for `frechet_mean(manifold_sphere(), x, ...)`
#' @param x an (n by q) array of data
#' @param tau step size
#' @param tol convergence tolerance
#' @param max.iter maximum number of iterations
#' @param init optional initial value
#' @param verbose if TRUE, print the loss at each iteration
#' @seealso [frechet_mean()]
#' @export
mean_on_sphere = function(x, tau = 0.1, tol = 1e-8, max.iter = 1000,
                          init = NULL, verbose = FALSE) {
  frechet_mean(manifold_sphere(), x, tau = tau, tol = tol, max.iter = max.iter,
              init = init, verbose = verbose)
}

#' Legacy alias for `parallel_transport(manifold_sphere(), x, y, V)`
#' @param x starting point
#' @param y end point
#' @param V tangent vector at x (a vector or k by d array)
#' @seealso [parallel_transport()]
#' @export
pt_sphere = function(x, y, V) parallel_transport(manifold_sphere(), x, y, V)
