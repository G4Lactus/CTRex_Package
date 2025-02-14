complex_cov <- function(X) {
  n <- nrow(X)
  X_mean <- scale(X, center = TRUE, scale = FALSE)
  cov_mat <- Conj(t(X)) %*% X

  return(cov_mat / (n - 1))
}
