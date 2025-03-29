#' @title
#' Generate Regularization Grid for Complex Ridge Regression
#'
#' @description
#' Internal function that constructs a logarithmically-spaced sequence of
#' regularization parameters (λ) for complex-valued ridge regression.
#' Creates λ values from ε·λ_max to λ_max for cross-validation and
#' regularization path computation.
#'
#' @param X Complex-valued predictor matrix (n x p).
#' @param y Complex-valued response vector (length n).
#' @param epsilon Ratio of smallest to largest λ (0 < ε < 1).
#' Determines λ_min = ε·λ_max.
#' @param L Number of λ values minus 1 (produces L+1 values).
#' Controls grid density.
#'
#' @return Numeric vector of (L+1) λ values sorted in increasing order.
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Generate complex data
#' X <- matrix(complex(real=rnorm(100), imaginary=rnorm(100)), ncol=10)
#' y <- complex(real=rnorm(10), imaginary=rnorm(10))
#'
#' # Create λ grid with 21 values (ε=0.01)
#' lambda_sequence <- lambda_grid(X, y, epsilon=0.01, L=20)
#'
#' # Visualize logarithmic spacing
#' plot(log(lambda_sequence), xlab="Index", ylab="log(λ)", type="b")
#'
#' }
#
lambda_grid <- function(X, y, epsilon, L) {

  n <- length(y)

  # Compute lambda_max
  lambdas <- apply(X, 2, function(col) { return( Mod(Conj(t(col)) %*% y) )})
  lambda_0 <- max(lambdas) / n

  # Generate lambda grid
  eps_power <- epsilon**(1 / L)
  lambda_values <- lambda_0 * eps_power**(seq.int(0, L))

  return(sort(lambda_values, decreasing = FALSE))
}
