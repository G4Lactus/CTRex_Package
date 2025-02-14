# Cross-validation to optimize lambda
lambda_cv <- function(X, y, epsilon, L, nfolds = 10) {

  n <- nrow(X)
  indices <- sample(rep(1:nfolds, length.out = n))
  lambdas <- lambda_grid(X, y, epsilon, L)

  cv_errors <- numeric(L)

  for (i in seq_along(lambdas)) {
    lambda <- lambdas[i]
    errors <- numeric(nfolds)

    for (k in 1:nfolds) {

      # Split into training and validation sets
      train_idx <- which(indices != k)
      val_idx <- which(indices == k)

      # Fit parameter vector beta
      beta <- complex_ridge(X[train_idx,], y[train_idx], lambda)

      # Predict on validation set
      errors[k] <- mean(Mod(y[val_idx] - X[val_idx,] %*% beta)^2)
    }

    cv_errors[i] <- mean(errors)
  }

  # Calculate the standard error of the cross-validation error
  se_errors <- sd(cv_errors) / sqrt(L)

  # Find the lambda with minimum MSE error
  best_lambda <- lambdas[which.min(cv_errors)]

  # Find lambda within 1se of the minimum error
  lambda_1se <- max(lambdas[cv_errors <= (min(cv_errors) + se_errors)])


  return(list(best_lambda = best_lambda,
              lambda_1se = lambda_1se,
              cv_errors = cv_errors,
              lambdas = lambdas))
}
