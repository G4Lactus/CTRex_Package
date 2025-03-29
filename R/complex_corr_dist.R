#' @title
#' Compute Complex Correlation Distance Matrix
#'
#' @description
#' Internal function that calculates a correlation-based distance matrix for
#'  complex-valued predictors.
#' Used in hierarchical clustering to group variables with similar correlation
#'  structures. The distance
#' between variables is defined as \(1 - |\rho|\), where \(\rho\) is the complex
#'  correlation coefficient.
#'
#' Key steps:
#' 1. Compute complex covariance matrix (without centering/scaling)
#' 2. Derive complex correlation matrix
#' 3. Convert correlations to distances using \(d = 1 - |\rho|\)
#'
#' @param X Complex-valued predictor matrix (n x p)
#'
#' @return stats::dist object containing pairwise distances between variables,
#'         suitable for hierarchical clustering
#'
#' @importFrom stats as.dist
#' @keywords internal
#'
#' @examples
#' \donttest{
#'
#' # Create complex-valued predictors
#' X <- matrix(complex(
#'   real = rnorm(30),
#'   imaginary = rnorm(30)),
#'   ncol = 5
#' )
#'
#' # Compute correlation distance matrix
#' dist_matrix <- complex_corr_dist(X)
#' print(as.matrix(dist_matrix)[1:3, 1:3])  # Show partial distance matrix
#'
#' }
#
complex_corr_dist <- function(X) {

  n <- nrow(X)

  # NOTE: scaling causes trouble with clustering a steering mtx.
  #X <- scale(X, center = TRUE, scale = FALSE)
  cov_mat <- Conj(t(X)) %*% X
  cov_mat <- cov_mat / (n - 1)

  # Compute variances (diagonal elements)
  vars <- diag(cov_mat)

  # Compute the correlation matrix
  corr_mat <- cov_mat / sqrt(outer(vars, vars, "*"))

  # Transform correlation matrix into distance matrix
  return(as.dist(1 - Mod(corr_mat)))

}
