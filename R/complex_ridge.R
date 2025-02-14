complex_ridge <- function(X, y, lambda) {
  p <- ncol(X)

  # Compute the solution
  XhX <- Conj(t(X)) %*% X
  Xhy <- Conj(t(X)) %*% y

  return(Matrix::solve(XhX + lambda * diag(p), Xhy))
}
