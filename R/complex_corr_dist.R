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
