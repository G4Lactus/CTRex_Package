# Compute lambda grid
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
