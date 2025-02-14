complex_hclust <- function(X, hc_method, coherence_max = 0.5) {

  # TODO: Data fidelity check

  # Perform hierarchical clustering
  hc_fit <- stats::hclust(complex_corr_dist(X), method = hc_method)
  clusters <- stats::cutree(hc_fit, h = 1 - coherence_max)
  max_clusters <- max(clusters)

  clusters <- data.frame(
    "Var" = names(clusters),
    "Cluster_Nr." = unname(clusters)
  )

  clusters <-
    stats::aggregate(clusters$"Var" ~ clusters$"Cluster_Nr.",
                     FUN = "c",
                     simplify = FALSE
    )


  cluster_sizes <- vector("numeric", length = max_clusters)
  for (j in seq(max_clusters)) {
    cluster_sizes[j] <- length(clusters$`clusters$Var`[[j]])
  }

  return(list(max_clusters = max_clusters,
              clusters = clusters,
              cluster_sizes = cluster_sizes))
}
