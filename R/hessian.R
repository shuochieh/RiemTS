#' Riemannian Hessian of the squared distance function
#'
#' Evaluates the action of the Riemannian Hessian of
#' \eqn{f(x) = 0.5 d^2(x, p)} at `x` on one or more tangent vectors `V`.
#'
#' @param mfd a manifold object
#' @param p target reference point
#' @param x base point where the Hessian is evaluated
#' @param V a (batch of) tangent vectors at `x`
#' @param ... additional arguments passed to methods
#'
#' @examples
#' mfd <- manifold_BWS()
#' x <- crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' p <- crossprod(matrix(rnorm(9), 3, 3)) + diag(1, 3)
#' V <- basis(mfd, x)
#' Hessian(mfd, p, x, V)
#'
#' @export
Hessian <- function(mfd, p, x, V, ...) UseMethod("Hessian")
