# RiemTS

Time series toolkit for data supported on Riemannian manifolds.

RiemTS provides a common interface for manifold operations — geodesic
distance, exponential/logarithm maps, parallel transport, Fréchet mean, and
the Riemannian Hessian — together with time series methods built on top of
them, such as the Riemannian factor model and change point detection.

Supported manifolds:

| Manifold | Constructor |
|---|---|
| Bures–Wasserstein (SPD matrices) | `manifold_BW()` |
| Log-Euclidean (SPD matrices) | `manifold_logE()` |
| Sphere | `manifold_sphere()` |
| Grassmannian (in progress) | `manifold_grassmann()` |

## Installation

```r
# install.packages("remotes")
remotes::install_github("shuochieh/RiemTS")
```

## Example

### Manifold operations

Every manifold shares the same generic functions, so code written for one
geometry works for another by swapping the manifold object.

```r
library(RiemTS)

mfd <- manifold_BW()

# two 3 x 3 SPD matrices
x <- crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
p <- crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)

geod(mfd, x, p)                # geodesic distance
v <- Log_mfd(mfd, x, p)        # tangent vector at x pointing to p
Exp_mfd(mfd, x, v)             # maps back to p

E <- basis(mfd, x)             # orthonormal basis of the tangent space at x
Riem_metric(mfd, x, v, v)      # squared length of v

# Riemannian Hessian of f(x) = 0.5 * d^2(x, p), applied to the basis vectors
Hessian(mfd, p, x, E)
```

### Fréchet mean

```r
x <- array(NA, dim = c(10, 3, 3))
for (i in 1:10) x[i, , ] <- crossprod(matrix(rnorm(9), 3, 3)) + diag(0.1, 3)

frechet_mean(mfd, x)                   # specialized algorithm
frechet_mean(mfd, x, method = "SGD")   # generic Riemannian SGD
```

### Riemannian factor model

`CaliforniaRain` holds the share of yearly precipitation falling in each
month (1895–2024). The square-root transform maps each year onto the sphere,
where we fit a two-factor model.

```r
data(CaliforniaRain)
x <- sqrt(CaliforniaRain)

fit <- rfm(manifold_sphere(), x, r = 2)
dim(fit$f_hat)   # 130 years x 2 factors
fit$E_val        # eigenvalues used to choose the number of factors
```

### Change point detection

`ZZS_WBS()` combines the test statistic of Zhang, Zhu, and Shao (2026) with
wild binary segmentation to locate multiple change points. It works with any
manifold object. Here the mean on the sphere jumps at time 100.

```r
set.seed(1)
mfd <- manifold_sphere()

x1 <- t(c(1, 0, 0) + t(matrix(rnorm(300, sd = 0.25), ncol = 3)))
x2 <- t(c(0, 0, 1) + t(matrix(rnorm(300, sd = 0.25), ncol = 3)))
x <- rbind(x1 / sqrt(rowSums(x1^2)), x2 / sqrt(rowSums(x2^2)))

U <- WBS_intervals(n = 200, M = 15)                # random sub-intervals
thres <- ZZS_WBS_thres_sim(n = 200, U = U)         # simulated 5% threshold
ZZS_WBS(mfd, x, U, thres = thres)                  # detected change points: 100
```

For a single change point, `ZZS_test()` returns the test statistic and the
location where it peaks.

## Roadmap

- [x] **Version 1**
  - [x] Double-check the Hessian in Log-Euclidean
- [ ] General-purpose Hessian approximation, and specialized formula for BW
- [ ] Manifold operations for the Grassmannian
- [ ] General-purpose MOSUM and CUSUM tests/detection
  - This step involves many tuning details
- [ ] **Version 2**

## References

Huang, S.-C., Chen, R., & Chen, Y. (2026). A Riemannian Factor Model for
Manifold-Valued Time Series. arXiv preprint arXiv:2607.28385.

Zhang, Y., Zhu, C., & Shao, X. (2026). Change-Point Detection for
Object-Valued Time Series. Journal of Business & Economic Statistics, 44(1),
255–269.
