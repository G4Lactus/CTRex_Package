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
