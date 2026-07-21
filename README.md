# manifoldstats

Statistics on matrix manifolds — Bures-Wasserstein SPD matrices, the
Grassmannian, and the sphere — through one shared interface, plus a
Riemannian factor model for Bures-Wasserstein data (`rfm_bws()`, from
`RFM_JBES.pdf`).

## Authors

- Shuo-chieh Huang
- Shen-Hsun Liao

Every manifold exposes the same five operations:

```r
geod(mfd, x, y)                    # geodesic distance
Exp_map(mfd, v, mu)                # exponential map (tangent space -> manifold)
Log_map(mfd, x, mu)                # logarithm map (manifold -> tangent space)
frechet_mean(mfd, x, ...)          # Riemannian (stochastic) gradient descent mean
parallel_transport(mfd, from, to, v)
```

`mfd` is a small object built by `manifold_bws()`, `manifold_grassmann(p)`, or
`manifold_sphere()` that tells these functions what a "point" looks like for
that manifold. Every function above works the same way whether `x`/`y`/`v`
is a single point or an `(n by ...)` batch of points.

## Installation

This is a local, unpublished package. From the project root:

```r
# for interactive development (no install step, re-run after edits)
pkgload::load_all(".")

# or a real install
install.packages(".", repos = NULL, type = "source")
library(manifoldstats)
```

Requires the `expm`, `deSolve`, and `maotai` packages (installed
automatically as dependencies).

## Quick start

### Pick a manifold

```r
mfd_bws    = manifold_bws()            # SPD matrices, Bures-Wasserstein metric
mfd_grass  = manifold_grassmann(p = 2) # Gr(d, 2): rank-2 projectors
mfd_sphere = manifold_sphere()         # unit sphere
```

### Distances, Exp/Log maps

Using the sphere as the simplest example — BWS and Grassmannian work
identically, just with matrix-valued points instead of vectors:

```r
mfd = manifold_sphere()
x = c(1, 0, 0)
y = c(0, 1, 0)

geod(mfd, x, y)                # pi/2

v = Log_map(mfd, x, y)         # tangent vector at y pointing towards x
Exp_map(mfd, v, y)             # maps back to (approximately) x
```

Batches work the same way — pass an `(n by q)` matrix (sphere) or an
`(n by p by p)` array (BWS/Grassmannian) instead of a single point, for
either `x`, `y`, or both:

```r
X = rbind(c(1, 0, 0), c(0, 0, 1), c(0, 1, 0))
geod(mfd, X, y)                 # length-3 vector of distances
Log_map(mfd, X, y)              # 3 x 3 matrix of tangent vectors at y
```

### Fréchet mean

```r
set.seed(1)
mfd = manifold_sphere()
X = matrix(rnorm(30 * 3), 30, 3)
X[, 1] = X[, 1] + 4             # cluster the points around c(1, 0, 0)
X = X / sqrt(rowSums(X^2))      # project onto the sphere

mu_hat = frechet_mean(mfd, X)

# stochastic gradient descent on a random subsample each iteration
frechet_mean(mfd, X, batch_size = 10)
```

`batch_size` is shared across all three manifolds: `NULL` means full-batch
gradient descent; a value in `(0, 1]` is a proportion of the data; a value
`> 1` is an absolute subsample count (both interpretations existed
separately in the pre-refactor code as `batch_por` and `batch_size`).

### Parallel transport

```r
mfd = manifold_sphere()
x = c(1, 0, 0)
y = c(0, 1, 0)
v = c(0, 0, 1)                  # tangent to the sphere at x
parallel_transport(mfd, x, y, v)
```

### Tangent-space bases and coordinates

BWS and the Grassmannian additionally support expressing a tangent vector
(a symmetric matrix) as a coordinate vector with respect to an orthonormal
basis, and reconstructing it:

```r
Sigma = diag(3) + 0.1
basis = tan_basis_bws(Sigma)$E            # orthonormal basis at Sigma

V = matrix(rnorm(9), 3, 3); V = (V + t(V)) / 2   # a symmetric tangent vector
z = coords_from_basis(basis, V)            # V's coordinates in that basis
tangent_from_coords(basis, z)              # back to V
```

The Grassmannian equivalent is `basis_gr(P, p)` with the same
`coords_from_basis()`/`tangent_from_coords()` pair.

## Bures-Wasserstein Riemannian Factor Model

`rfm_bws()` fits the Riemannian factor model of `RFM_JBES.pdf` to a sample
of SPD matrices: it estimates the Fréchet mean, builds a tangent-space
basis there, and fits a linear factor model (`LYB_fm()`, Lam-Yao-Bathia) on
the log-mapped, coordinate-expressed data.

```r
set.seed(1)
n = 60; p = 4
X = array(NA, dim = c(n, p, p))
for (i in 1:n) {
  A = matrix(rnorm(p * p), p, p)
  X[i, , ] = A %*% t(A) + diag(p)
}

model = rfm_bws(X, r = 2, h = 6)
model$mu_hat            # estimated Frechet mean
dim(model$f_hat)        # estimated factor process, n x r
model$A                 # loading matrix

# evaluate fraction of variance explained by 1..r factors on held-out data
Frac_Var_bws(X[1:10, , ], model, Euclidean_mean = model$mu_hat)
```

`Frac_Var_LYB()` evaluates a plain (non-Riemannian) linear factor model
fit with `LYB_fm()` directly on vectorized SPD matrices, for comparison.
`dyn_RFM()`/`dyn_LFM()` build one-step-ahead forecasts using a VAR(1) on
the estimated factors.

> **Note:** `dyn_RFM()` calls `main_BWS(...)`, which isn't defined
> anywhere in this package — that's a pre-existing gap from before this
> package existed, not something introduced here.

## Package layout

| File | Contents |
|---|---|
| `R/manifold.R` | manifold constructors, generic dispatch machinery, `geod`/`Exp_map`/`Log_map`, `coords_from_basis`/`tangent_from_coords` |
| `R/bws.R` | Bures-Wasserstein core math + S3 methods, `tan_basis_bws`, `log_vec_construct` |
| `R/grassmann.R` | Grassmannian core math + S3 methods, `basis_gr` |
| `R/sphere.R` | sphere core math + S3 methods, `basis_sphere` |
| `R/frechet_mean.R` | generic Fréchet mean (Riemannian gradient descent) |
| `R/parallel_transport.R` | generic `parallel_transport()` dispatcher |
| `R/compat.R` | old per-manifold function names (`geod_BWS`, `mean_on_gr`, `Exp_sphere`, ...), kept as one-line wrappers |
| `R/bws_factor_model.R` | `rfm_bws`, `LYB_fm`, `Frac_Var_*`, `dyn_RFM`/`dyn_LFM`, and related BWS-specific modeling code |

Every exported function has full documentation with runnable examples —
see `?geod`, `?frechet_mean`, `?rfm_bws`, etc., or browse the generated
`man/*.Rd` files. For the history and rationale of the refactor that
produced this structure, see `CHANGES.html`.

## Tests

```sh
Rscript tests/smoke_test.R
```

A dependency-free script (no `testthat`) checking Exp/Log round-trips,
`geod(x,x) == 0`, that `frechet_mean()` numerically matches the original
pre-refactor per-manifold implementations, and `parallel_transport()` for
all three manifolds.
