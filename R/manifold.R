################################################################################
#### Generic manifold interface: constructors + shared array-dispatch helpers #
################################################################################
#
# The three manifolds implemented in this package (Bures-Wasserstein SPD
# matrices, the Grassmannian, the sphere) all expose the same primitives
# (geod, Exp_map, Log_map, frechet_mean, parallel_transport) but represent a
# "single point" differently: BWS/Grassmannian points are p-by-p matrices
# (point_ndim = 2), sphere points are length-q vectors (point_ndim = 1). A
# "batch of n points" always has one extra leading dimension. Everything
# below is generic in that one fact, so the single-vs-batch dispatch logic
# that used to be copy-pasted in every geod_*/Exp_*/Log_* wrapper only needs
# to be written once.

#' @importFrom maotai lyapunov
#' @importFrom expm sqrtm expm
#' @importFrom deSolve ode
NULL

#' Bures-Wasserstein manifold of SPD matrices
#'
#' Construct a manifold object to pass as the first argument of [geod()],
#' [Exp_map()], [Log_map()], [frechet_mean()], and [parallel_transport()].
#' Points on this manifold are symmetric positive-definite (p by p) matrices.
#'
#' @return an object of class `c("manifold_bws", "manifold")`
#' @examples
#' mfd = manifold_bws()
#'
#' # two random 3x3 SPD matrices
#' A = matrix(rnorm(9), 3, 3); A = A %*% t(A) + diag(3)
#' B = matrix(rnorm(9), 3, 3); B = B %*% t(B) + diag(3)
#'
#' geod(mfd, A, B)                 # geodesic (Bures-Wasserstein) distance
#' v = Log_map(mfd, B, A)          # tangent vector at A pointing towards B
#' Exp_map(mfd, v, A)              # maps back to (approximately) B
#' @export
manifold_bws = function() {
  structure(list(point_ndim = 2), class = c("manifold_bws", "manifold"))
}

#' Grassmannian Gr(d, p)
#'
#' Construct a manifold object to pass as the first argument of [geod()],
#' [Exp_map()], [Log_map()], [frechet_mean()], and [parallel_transport()].
#' Points on this manifold are rank-p orthogonal projectors (d by d matrices).
#'
#' @param p rank of the projectors; if NULL it is inferred per-call from the data
#' @return an object of class `c("manifold_grassmann", "manifold")`
#' @examples
#' d = 5; p = 2
#' mfd = manifold_grassmann(p)
#'
#' # two random rank-p projectors
#' rand_projector = function(d, p) {
#'   U = qr.Q(qr(matrix(rnorm(d * d), d, d)))[, 1:p, drop = FALSE]
#'   U %*% t(U)
#' }
#' A = rand_projector(d, p)
#' B = rand_projector(d, p)
#'
#' geod(mfd, A, B)
#' @export
manifold_grassmann = function(p = NULL) {
  structure(list(point_ndim = 2, p = p), class = c("manifold_grassmann", "manifold"))
}

#' Sphere
#'
#' Construct a manifold object to pass as the first argument of [geod()],
#' [Exp_map()], [Log_map()], [frechet_mean()], and [parallel_transport()].
#' Points on this manifold are unit-norm vectors.
#'
#' @return an object of class `c("manifold_sphere", "manifold")`
#' @examples
#' mfd = manifold_sphere()
#' x = c(1, 0, 0)
#' y = c(0, 1, 0)
#'
#' geod(mfd, x, y)                 # pi / 2, a quarter turn apart
#'
#' # batches of points work the same way, as an (n by q) matrix
#' X = rbind(c(1, 0, 0), c(0, 0, 1))
#' geod(mfd, X, y)
#' @export
manifold_sphere = function() {
  structure(list(point_ndim = 1), class = c("manifold_sphere", "manifold"))
}

is_batch = function(mfd, x) {
  if (mfd$point_ndim == 1) {
    return(is.matrix(x))
  } else {
    return(length(dim(x)) == mfd$point_ndim + 1)
  }
}

get_point = function(mfd, x, i) {
  if (mfd$point_ndim == 1) {
    return(x[i, ])
  } else {
    return(x[i, , ])
  }
}

n_points = function(mfd, x) {
  if (mfd$point_ndim == 1) {
    return(nrow(x))
  } else {
    return(dim(x)[1])
  }
}

subset_points = function(mfd, x, idx) {
  if (mfd$point_ndim == 1) {
    return(x[idx, , drop = FALSE])
  } else {
    return(x[idx, , , drop = FALSE])
  }
}

# average a batch of tangent vectors/matrices produced by Log_map()
mean_tangent = function(mfd, v) {
  if (!is_batch(mfd, v)) {
    return(v)
  }
  if (mfd$point_ndim == 1) {
    return(colMeans(v))
  } else {
    return(apply(v, seq_len(mfd$point_ndim) + 1, mean))
  }
}

#' Map a single-point function over a batch of points
#'
#' Replaces the repeated "is x a single point or an (n by ...) array of
#' points?" dispatch that used to live in every Exp_*/Log_* wrapper.
#'
#' @param mfd a manifold object
#' @param core_fn a function taking a single point and returning either a
#'   point (same family as the input) or a coordinate vector
#' @param x a single point or a batch of points
map_over_points = function(mfd, core_fn, x) {
  if (!is_batch(mfd, x)) {
    return(core_fn(x))
  }

  n = n_points(mfd, x)
  first = core_fn(get_point(mfd, x, 1))

  if (is.matrix(first)) {
    res = array(NA, dim = c(n, dim(first)))
    res[1, , ] = first
    if (n > 1) {
      for (i in 2:n) res[i, , ] = core_fn(get_point(mfd, x, i))
    }
  } else {
    res = matrix(NA, n, length(first))
    res[1, ] = first
    if (n > 1) {
      for (i in 2:n) res[i, ] = core_fn(get_point(mfd, x, i))
    }
  }

  return(res)
}

#' Map a pairwise (point, point) -> scalar function over batches
#'
#' Replaces the repeated (single,single)/(array,single)/(single,array)/
#' (array,array) branching that used to live in every geod_* wrapper.
#'
#' @param mfd a manifold object
#' @param core_fn a function taking two single points and returning a scalar
#' @param x,y single points or batches of points
map_over_point_pairs = function(mfd, core_fn, x, y) {
  x_batch = is_batch(mfd, x)
  y_batch = is_batch(mfd, y)

  if (!x_batch && !y_batch) {
    return(core_fn(x, y))
  } else if (x_batch && !y_batch) {
    n = n_points(mfd, x)
    res = rep(NA, n)
    for (i in 1:n) res[i] = core_fn(get_point(mfd, x, i), y)
  } else if (!x_batch && y_batch) {
    n = n_points(mfd, y)
    res = rep(NA, n)
    for (i in 1:n) res[i] = core_fn(x, get_point(mfd, y, i))
  } else {
    nx = n_points(mfd, x)
    ny = n_points(mfd, y)
    if (nx != ny) {
      stop("map_over_point_pairs: x and y have different sample sizes")
    }
    res = rep(NA, nx)
    for (i in 1:nx) res[i] = core_fn(get_point(mfd, x, i), get_point(mfd, y, i))
  }

  return(res)
}

geod_core = function(mfd, x, y, ...) UseMethod("geod_core")
Exp_core = function(mfd, v, mu, ...) UseMethod("Exp_core")
Log_core = function(mfd, x, mu, ...) UseMethod("Log_core")

#' Geodesic distance on a manifold
#'
#' Works identically for [manifold_bws()], [manifold_grassmann()], and
#' [manifold_sphere()] objects, and for any combination of single points and
#' batches of points (an (n by ...) array/matrix) in `x` and `y`.
#'
#' @param mfd a manifold object (manifold_bws(), manifold_grassmann(p), manifold_sphere())
#' @param x,y single points or (n by ...) batches of points
#' @return a scalar (single point vs. single point), or a length-n vector if
#'   either `x` or `y` (or both, paired) is a batch of n points
#' @examples
#' mfd = manifold_sphere()
#' x = c(1, 0, 0)
#' y = c(0, 1, 0)
#' geod(mfd, x, y)
#'
#' # one point vs. a batch of points
#' X = rbind(c(1, 0, 0), c(0, 0, 1), c(0, 1, 0))
#' geod(mfd, X, y)
#' @export
geod = function(mfd, x, y) {
  map_over_point_pairs(mfd, function(a, b) geod_core(mfd, a, b), x, y)
}

#' Exponential map on a manifold
#'
#' Maps a tangent vector `v` at base point `mu` back onto the manifold.
#' Inverse of [Log_map()].
#'
#' @param mfd a manifold object
#' @param v a single tangent vector or a batch of tangent vectors, at mu
#' @param mu base point
#' @return a point, or a batch of points, in the same shape as `v`
#' @examples
#' mfd = manifold_sphere()
#' mu = c(1, 0, 0)
#' x = c(0, 1, 0)
#' v = Log_map(mfd, x, mu)   # tangent vector at mu pointing towards x
#' Exp_map(mfd, v, mu)       # back to (approximately) x
#' @export
Exp_map = function(mfd, v, mu) {
  map_over_points(mfd, function(vi) Exp_core(mfd, vi, mu), v)
}

#' Logarithm map on a manifold
#'
#' Maps a point `x` to the tangent space at base point `mu`. Inverse of
#' [Exp_map()].
#'
#' @param mfd a manifold object
#' @param x a single point or a batch of points
#' @param mu base point
#' @return a tangent vector, or a batch of tangent vectors, at `mu`
#' @examples
#' mfd = manifold_sphere()
#' mu = c(1, 0, 0)
#' x = c(0, 1, 0)
#' Log_map(mfd, x, mu)
#' @export
Log_map = function(mfd, x, mu) {
  map_over_points(mfd, function(xi) Log_core(mfd, xi, mu), x)
}

#' Coordinates of tangent matrix/matrices w.r.t. a basis
#'
#' Unifies the (mathematically identical) coord_gr()/tangent_in_E() from the
#' old Grassmannian/BWS utilities: for each basis element E_i, computes
#' 0.5 * trace(E_i %*% x).
#'
#' @param basis an (m by p by p) array of basis matrices
#' @param x a (p by p) matrix or an (n by p by p) array
#' @return a length-m coordinate vector, or an (n by m) matrix if `x` is a batch
#' @examples
#' Sigma = diag(3) + 0.1
#' basis = tan_basis_bws(Sigma)$E   # an orthonormal basis for the tangent space at Sigma
#'
#' V = matrix(rnorm(9), 3, 3); V = (V + t(V)) / 2   # a random symmetric tangent vector
#' z = coords_from_basis(basis, V)                  # V's coordinates in that basis
#' V_reconstructed = tangent_from_coords(basis, z)
#' max(abs(V - V_reconstructed))                    # ~ 0
#' @export
coords_from_basis = function(basis, x) {
  m = dim(basis)[1]

  if (length(dim(x)) == 3) {
    n = dim(x)[1]
    res = matrix(NA, n, m)
    for (j in 1:n) {
      for (i in 1:m) {
        res[j, i] = 0.5 * sum(diag(basis[i, , ] %*% x[j, , ]))
      }
    }
  } else {
    res = rep(NA, m)
    for (i in 1:m) {
      res[i] = 0.5 * sum(diag(basis[i, , ] %*% x))
    }
  }

  return(res)
}

#' Reconstruct tangent matrix/matrices from coordinates w.r.t. a basis
#'
#' Unifies the (mathematically identical) coord2tan()/log_to_tangent() from
#' the old Grassmannian/BWS utilities.
#'
#' @param basis an (m by p by p) array of basis matrices
#' @param coords a length-m coefficient vector, or an (n by m) matrix
#' @return a (p by p) matrix, or an (n by p by p) array if `coords` is a matrix
#' @examples
#' Sigma = diag(3) + 0.1
#' basis = tan_basis_bws(Sigma)$E
#' z = c(1, 0, 0, 0, 0, 0)          # coordinates of the first basis element
#' tangent_from_coords(basis, z)    # reconstructs basis[1, , ]
#' @export
tangent_from_coords = function(basis, coords) {
  p = dim(basis)[2]

  if (is.matrix(coords)) {
    n = nrow(coords)
    m = ncol(coords)
    res = array(0, dim = c(n, p, p))
    for (j in 1:n) {
      for (i in 1:m) {
        res[j, , ] = res[j, , ] + coords[j, i] * basis[i, , ]
      }
    }
  } else {
    res = matrix(0, p, p)
    for (i in seq_along(coords)) {
      res = res + coords[i] * basis[i, , ]
    }
  }

  return(res)
}
