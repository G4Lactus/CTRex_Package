#' True positive proportion (TPP)
#'
#' Computes the TPP based on the estimated and the true regression coefficient vectors.
#'
#' @param beta_hat Estimated regression coefficient vector.
#' @param beta True regression coefficient vector.
#' @param eps Numerical zero.
#'
#' @return True positive proportion (TPP).
#'
#' @export
#'
#' @examples
#' set.seed(2)
# ----------------------------------------------------------------------
TPP <- function(beta_hat,
                beta,
                eps = .Machine$double.eps) {
  # Remove all dimension attributes of length one
  beta_hat <- drop(beta_hat)
  beta <- drop(beta)

  # Error control
  if (!is.vector(beta_hat)) {
    stop("'beta_hat' must be a vector.")
  }

  if (!is.numeric(beta_hat)) {
    stop("'beta_hat' only allows numerical values.")
  }

  if (anyNA(beta_hat)) {
    stop("'beta_hat' contains NAs. Please remove or impute them before proceeding.")
  }

  if (!is.vector(drop(beta))) {
    stop("'beta' must be a vector.")
  }

  if (!is.numeric(beta)) {
    stop("'beta' only allows numerical values.")
  }

  if (anyNA(beta)) {
    stop("'beta' contains NAs. Please remove or impute them before proceeding.")
  }

  if (length(beta_hat) != length(beta)) {
    stop("Length of beta_hat does not match length of beta.")
  }

  # Compute TPP
  num_actives <- sum(abs(beta) > eps)
  num_true_positives <- sum(abs(beta) > eps & abs(beta_hat) > eps)

  if (num_actives == 0) {
    tpp <- 0
  } else {
    tpp <- (num_true_positives / num_actives)
  }
  return(tpp)
}
