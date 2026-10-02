# Generic manifold interface: constructors + shared array-dispatch helpers 

#' @importFrom maotai lyapunov
#' @importFrom expm sqrtm expm
#' @importFrom deSolve ode
NULL

#' Sphere
#'
#' Construct a manifold object. Points on this manifold are unit-norm vectors.
#'
#' @return an object of class `c("manifold_sphere", "manifold")`
#' @export
manifold_sphere = function() {
  structure(list(point_ndim = 1), class = c("manifold_sphere", "manifold"))
}

#' Bures-Wasserstein manifold of SPD matrices
#'
#' Construct a manifold object. Points on this manifold are symmetric 
#' positive-definite (m by m) matrices.
#'
#' @return an object of class `c("manifold_BWS", "manifold")`
#' @export
manifold_bws = function() {
  structure(list(point_ndim = 2), class = c("manifold_BWS", "manifold"))
}

#' Log-Euclidean manifold of SPD matrices
#'
#' Construct a manifold object. Points on this manifold are symmetric 
#' positive-definite (m by m) matrices.
#'
#' @return an object of class `c("manifold_logEuclidean", "manifold")`
#' @export
manifold_logE = function() {
  structure(list(point_ndim = 2), class = c("manifold_logE", "manifold"))
}

#' Grassmannian Gr(r, d)
#'
#' Construct a manifold object. Points on this manifold are rank-r orthogonal 
#' projectors (d by d matrices).
#'
#' @param r rank of the projectors; if NULL it is inferred per-call from the data
#' @return an object of class `c("manifold_grassmann", "manifold")`
#' @export
manifold_grassmann = function (r = NULL) {
  structure(list(point_ndim = 2, r = r), class = c("manifold_grassmann", "manifold"))
}

n_points = function (mfd, x) {
  if (mfd$point_ndim == 1) {
    return (nrow(x))
  } else {
    return (dim(x)[1])
  }
}

get_point = function (mfd, x, i) {
  if (mfd$point_ndim == 1) {
    return (x[i,])
  } else {
    return (x[i,,])
  }
}

subset_points = function (mfd, x, idx) {
  if (mfd$point_ndim == 1) {
    return(x[idx, , drop = FALSE])
  } else {
    return(x[idx, , , drop = FALSE])
  }
}

#' Geodesic distance on a manifold
#' 
#' Works identically for [manifold_sphere()], [manifold_BWS()],  
#' [manifold_logEuclidean()], and [manifold_grassmann()] objects, and for any 
#' combination of single points and batches of points (an (n by ...) array/matrix) 
#' in `x` and `y`.
#' 
#' @param mfd a manifold object
#' @param x single point or (n by ...) batches of points
#' @param y single point or (n by ...) batches of points
#' @return a scalar (single point vs. single point), or a length-\eqn{n} vector if
#'   either `x` or `y` (or both, paired) is a batch of \eqn{n} points
#'   
#' @examples
#' mfd = manifold_sphere()
#' y = c(0, 1, 0)
#' X = rbind(c(1, 0, 0), c(0, 0, 1), c(0, 1, 0))
#' geod(mfd, X, y)
#' 
#' @export
geod = function (mfd, x, y, ...) UseMethod("geod")

#' Exponential map on a manifold
#'
#' Maps a tangent vector `v` at base point `p` to the manifold.
#'
#' @param mfd a manifold object
#' @param p base point
#' @param v a single tangent vector or a batch of tangent vectors, at `p`
#' @return a point, or a batch of points, in the same shape as `v`
#' 
#' @examples
#' mfd = manifold_sphere()
#' p = c(1, 0, 0)
#' v = c(0, 0.5, 0)
#' Exp_mfd(mfd, p, v)       
#' 
#' @export
Exp_mfd = function (mfd, p, v, ...) UseMethod("Exp_mfd")

#' Logarithm map on a manifold
#'
#' Maps a point `q` to the tangent space at base point `p`. 
#'
#' @param mfd a manifold object
#' @param p base point
#' @param q a single point or a batch of points
#' @return a tangent vector, or a batch of tangent vectors, at `p`
#' 
#' @examples
#' mfd = manifold_sphere()
#' p = c(1, 0, 0)
#' q = c(0, 1, 0)
#' Log_mfd(mfd, p, q)
#' 
#' @export
Log_mfd = function (mfd, p, q, ...) UseMethod("Log_mfd")

#' Parallel transport on a manifold
#' 
#' Parallel-transport a (batch of) tangent vectors `v` 
#' 
#' @param mfd a manifold object
#' @param from starting point on the manifold
#' @param to end point of the transport
#' @param v a (batch of) tangent vectors at `from`
#' @return a (batch of) tangent vectors at `to`
#' 
#' @examples 
#' mfd = manifold_sphere()
#' from = c(1, 0, 0)
#' to = c(0, 1, 0)
#' v = c(0, 0.5, 0)
#' ptransport(mfd, from, to, v)
#' 
#' @export
ptransport = function (mfd, from, to, v, ...) UseMethod("ptransport")

#' Evaluate the Riemannian metric between two (batches of) tangent vectors at a
#' point `p`
#' 
#' @param mfd a manifold object
#' @param p base point
#' @param v a (batch of) tangent vectors
#' @param w a (batch of) tangent vectors
#' 
#' @examples 
#' mfd = manifold_sphere()
#' p = c(1, 0, 0)
#' v = c(0, 1, 0)
#' w = c(0, 0, 1)
#' Riem_metric(mfd, p, v, w)
#' 
#' @export
Riem_metric = function (mfd, p, v, w, ...) UseMethod("Riem_metric")




