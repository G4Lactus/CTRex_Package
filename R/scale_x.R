scale_x <- function(x,
                    num_rows = nrow(x),
                    num_cols = ncol(x)) {

  # center values
  x <- scale(x, center = TRUE, scale = FALSE)

  # determine column standard deviations
  x_cstds <- apply(x, 2, function(xs) {
    return(sqrt(sum(Mod(xs) ** 2)  / (num_rows - 1)))
  })

  # normalize
  for (colX in seq_len(num_cols)) {
    x[, colX] <- x[, colX] / x_cstds[colX]
  }
  return(x)
}