#' @title
#' Compute Complex-Valued Covariance Matrix
#'
#' @description
#' Internal function that calculates the sample covariance matrix for complex-valued data.
#' Centers variables by column means and applies Bessel's correction (n-1 denominator).
#' Core component for correlation/distance calculations in CT-Rex workflows.
#'
#' @param X Complex-valued predictor matrix (n x p)
#'
#' @return Complex covariance matrix (p x p) where:
#' \itemize{
#'   \item Diagonal elements are real-valued variances
#'   \item Off-diagonal elements are complex covariances
#' }
#'
#' @importFrom stats scale
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Create complex data matrix
#' X <- matrix(complex(
#'   real = rnorm(20),
#'   imaginary = rnorm(20)),
#'   ncol = 4
#' )
#'
#' # Compute covariance matrix
#' cov_mat <- complex_cov(X)
#' print(cov_mat[1:2, 1:2])  # Show top-left submatrix
#' }
#'
complex_cov <- function(X) {

  n <- nrow(X)
  X_mean <- scale(X, center = TRUE, scale = FALSE)
  cov_mat <- Conj(t(X)) %*% X

  return(cov_mat / (n - 1))
}
