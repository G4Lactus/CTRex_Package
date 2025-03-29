#' @title
#' Complex-Valued Ridge Regularization Tuning via Cross-Validation
#'
#' @description
#' Internal function performing k-fold cross-validation to select optimal
#' L2 regularization parameter (λ) for complex ridge regression. Implements
#' both minimum MSE and 1-standard-error rule for λ selection.
#'
#' @param X Complex-valued predictor matrix (n x p).
#' @param y Complex-valued response vector (length n).
#' @param epsilon Ratio of smallest to largest λ in grid (0 < ε < 1).
#' @param L Number of λ values to evaluate in regularization path.
#' @param nfolds Number of cross-validation folds. Default: 10.
#'
#' @return List containing: \itemize{
#'   \item best_lambda: λ achieving minimum cross-validated MSE
#'   \item lambda_1se: Largest λ within 1 standard error of minimum MSE
#'   \item cv_errors: Average MSE for each λ (length L)
#'   \item lambdas: Tested λ sequence (length L)
#' }
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Generate complex data
#' X <- matrix(complex(real=rnorm(100), imaginary=rnorm(100)), ncol=10)
#' y <- complex(real=rnorm(10), imaginary=rnorm(10))
#'
#' # Find optimal lambda with 20 values and 5 folds
#' cv_result <- lambda_cv(X, y, epsilon=0.01, L=20, nfolds=5)
#'
#' # Plot CV curve
#' plot(log(cv_result$lambdas), cv_result$cv_errors, type="b")
#' abline(v=log(cv_result$best_lambda), col="red")
#'
#' }
#
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
