#' @title
#' Complex-Valued Ridge Regression Solution
#'
#' @description
#' Internal function computing ridge regression estimates for complex-valued data.
#' Solves the regularized linear system using L2 penalty. Core component for
#' regularization parameter estimation in CT-Rex's group selection methods.
#'
#' @param X Complex-valued predictor matrix (n x p).
#' @param y Complex-valued response vector (length n).
#' @param lambda Positive regularization parameter controlling shrinkage strength.
#'
#' @return Complex-valued coefficient vector of length p.
#'
#' @importFrom Matrix solve
#' @importFrom stats Conj
#' @keywords internal
#'
#' @examples
#' \donttest{
#'
#' # Generate complex data
#' X <- matrix(complex(real=rnorm(20), imaginary=rnorm(20)), ncol=2)
#' y <- complex(real=rnorm(10), imaginary=rnorm(10))
#'
#' # Compute ridge solution
#' beta <- complex_ridge(X, y, lambda = 0.1)
#' print(Re(beta)) # Real components
#' print(Im(beta)) # Imaginary components
#'
#' }
#
complex_ridge <- function(X, y, lambda) {
  p <- ncol(X)

  # Compute the solution
  XhX <- Conj(t(X)) %*% X
  Xhy <- Conj(t(X)) %*% y

  return(Matrix::solve(XhX + lambda * diag(p), Xhy))
}
