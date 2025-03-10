update_cmplx_gvs_dummies <- function(cenv = environment()) {

  X_p_sub_dummy <- matrix(NA, nrow = cenv$n, ncol = cenv$p)

  p <- cenv$p
  for (w in seq(cenv$w_max)) {
    for (z in seq(cenv$cluster_res$max_clusters)) {

      sub_X <- cenv$X[, cenv$cluster_res$clusters$`clusters$Var`[[z]], drop = FALSE]
      sigma_sub_X <- complex_cov(sub_X)
      idx <- cumsum(cenv$cluster_res$cluster_sizes)
      mu <- rep(0, times = cenv$cluster_res$cluster_sizes[[z]])

      if (z == 1) {
        X_p_sub_dummy[, seq(idx[z])] <- cmvnorm::rcmvnorm(
                                          cenv$n,
                                          mean = mu,
                                          sigma = sigma_sub_X)

      } else {
        X_p_sub_dummy[, seq(idx[z - 1] + 1, idx[z])] <- cmvnorm::rcmvnorm(
                                                          cenv$n,
                                                          mean = mu,
                                                          sigma = sigma_sub_X)

      }
    }
    cenv$X_dummy[, seq(w * p + 1, (w + 1) * p)] <- X_p_sub_dummy
  }


}
