update_cmplx_gvs_dummies <- function(cenv = environment()) {

  n <- cenv$n
  p <- cenv$p
  w_max <- cenv$num_dummies / p
  X_p_sub_dummy_blck <- matrix(NA, nrow = n, ncol = p)

  cluster_sizes <- cenv$cluster_res$cluster_sizes
  cluster_var <- cenv$cluster_res$clusters$`clusters$Var`
  max_clusters <- cenv$cluster_res$max_clusters
  idx <- cumsum(cluster_sizes)
  X_dummy <- cenv$X_dummy

  for (w in seq_len(w_max)) {
    for (z in seq(max_clusters)) {

      sigma_sub_X <- complex_cov(cenv$X[, cluster_var[[z]], drop = FALSE])
      mu <- rep(0, cluster_sizes[z])

      if (z == 1) {
        X_p_sub_dummy_blck[, seq(idx[z])] <- cmvnorm::rcmvnorm(
                                          n = n,
                                          mean = mu,
                                          sigma = sigma_sub_X)

      } else {
        X_p_sub_dummy_blck[, seq(idx[z - 1] + 1, idx[z])] <- cmvnorm::rcmvnorm(
                                                          n = n,
                                                          mean = mu,
                                                          sigma = sigma_sub_X)

      }
    }
    X_dummy[, seq.int(w * p + 1, (w + 1) * p)] <- X_p_sub_dummy_blck
  }
  cenv$X_dummy <- X_dummy
}
