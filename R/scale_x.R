#' @title
#' Complex Matrix Centering and Scaling
#'
#' @description
#' Internal function that centers and scales complex-valued matrices. Performs:
#' 1. Column-wise mean centering
#' 2. Scaling by magnitude-based standard deviations
#' 3. Protection against near-constant columns via numerical stability threshold
#'
#' @param x Complex-valued matrix (n x p)
#' @param num_rows Number of observations. Default: `nrow(x)`
#' @param num_cols Number of variables. Default: `ncol(x)`
#'
#' @return Centered and scaled complex matrix with: \itemize{
#'   \item Column means = 0+0i
#'   \item Column norms scaled to sqrt(n) unless near-constant
#' }
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Create complex matrix with one constant column
#' X <- matrix(complex(
#'     real = c(rnorm(90), rep(0,10)),
#'     imaginary = c(rnorm(90), rep(0,10))),
#'   ncol = 10
#' )
#'
#' # Scale matrix
#' X_scaled <- scale_x(X)
#'
#' # Verify centering
#' print(apply(Re(X_scaled), 2, mean))  # Near-zero real means
#' print(apply(Im(X_scaled), 2, mean))  # Near-zero imaginary means
#' }
#
scale_x <- function(x,
                    num_rows = nrow(x),
                    num_cols = ncol(x)) {

  # center values
  x <- scale(x, center = TRUE, scale = FALSE)

  if (num_rows < 2) {
    denum <- num_rows
  } else {
    denum <- num_rows - 1
  }
  sqrt_n <- sqrt(num_rows)

  # determine column standard deviations
  x_cstds <- apply(x, 2, function(xs) {
    return(sqrt(sum(Mod(xs) ** 2)  / denum))
  })

  # normalize
  tol <- .Machine$double.eps
  for (colX in seq_len(num_cols)) {
    col_std <- x_cstds[colX]
    if (col_std / sqrt_n < tol) {
      x[, colX] <- tol * sqrt_n
    } else {
      x[, colX] <- x[, colX] / col_std
    }
  }

  return(x)
}
