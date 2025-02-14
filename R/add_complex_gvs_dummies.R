add_complex_gvs_dummies <- function(X,
                                    num_dummies,
                                    gvs_type,
                                    hc_method,
                                    coherence_max = 0.5) {

  # TODO: data fidelity check (see below - open)
  n <- nrow(X)
  p <- ncol(X)
  colnames(X) <- paste0("V", seq(p))

  # Perform clustering
  cluster_res <- complex_hclust(X, hc_method = hc_method, coherence_max = coherence_max)

  # Generate dummy predictors and append them to original X
  w_max <- num_dummies / p
  X_dummy <- matrix(NA, nrow = n, ncol = p + num_dummies)
  X_dummy[, seq(p)] <- X

  if (gvs_type == "cIEN") {
    # Binary cluster identity vectors for T-Rex+GVS+IEN
    # ----------------------------------------------------
    IEN_cl_ind <- lapply(cluster_res$clusters$`clusters$Var`, FUN = function(x) as.numeric(sub(".", "", x)))
    IEN_cl_id_vectors <- t(sapply(seq_along(IEN_cl_ind), FUN = function(x) {
      cl_id_vec <- rep(FALSE, times = ncol(X))
      cl_id_vec[IEN_cl_ind[[x]]] <- TRUE
      cl_id_vec}))
  }

  # Update complex gvs dummies
  update_cmplx_gvs_dummies(cenv = environment())

  # Output
  # ------------------------------------
  if (gvs_type == "cEN") { # ------------------------------
    return(list(
      X_dummy = X_dummy,
      max_clusters = cluster_res$max_clusters,
      cluster_sizes = cluster_res$cluster_sizes
    ))

  } else if (gvs_type == "cIEN") { # ----------------------
    return(
      list(X_dummy = X_dummy,
           max_clusters = cluster_res$max_clusters,
           cluster_sizes = cluster_res$cluster_sizes,
           IEN_cl_id_vectors = IEN_cl_id_vectors
      ))

  } else {
    stop("`gvs_type` not supported.")
  }

}


# TODO: data fidelity check

