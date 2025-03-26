complex_hclust <- function(X,
                           hc_method,
                           tree_cut_type,
                           coherence_max = 0.5,
                           minModuleSize = 10) {

  # TODO: Input data fidelity check

  # Perform hierarchical clustering
  hc_fit <- stats::hclust(complex_corr_dist(X), method = hc_method)

  if (tree_cut_type == "fixed") {
    clusters <- stats::cutree(hc_fit, h = 1 - coherence_max)
    max_clusters <- max(clusters)

  } else if (tree_cut_type == "dynamic") {
    if (is.null(minModuleSize)) {
      # Perform model order determination
      minModuleSizes <- seq(5, 50, by = 5)
      optimal_minModuleSize <- determine_optimal_minModuleSize(hc_fit, minModuleSizes)
      minModuleSize <- optimal_minModuleSize

    }
    clusters <- dynamicTreeCut::cutreeDynamicTree(
                        hc_fit,
                        deepSplit = FALSE,
                        minModuleSize = minModuleSize
                      )
    clusters <- clusters + 1
    names(clusters) <- paste0("V", 1:length(clusters))
    max_clusters <- max(clusters)

  } else {
    stop("`tree_cut_type` not supported.")
  }


  clusters <- data.frame(
    "Var" = names(clusters),
    "Cluster_Nr." = unname(clusters)
  )

  clusters <- stats::aggregate(clusters$"Var" ~ clusters$"Cluster_Nr.",
                     FUN = "c",
                     simplify = FALSE)

  cluster_sizes <- vector("numeric", length = max_clusters)
  for (j in seq(max_clusters)) {
    cluster_sizes[j] <- length(clusters$`clusters$Var`[[j]])
  }

  return(
    list(max_clusters = max_clusters,
         clusters = clusters,
         cluster_sizes = cluster_sizes
         )
  )

}



determine_optimal_minModuleSize <- function(hc_fit, minModuleSizes) {
  # Initialize vectors to store metrics
  nClusters <- numeric(length(minModuleSizes))

  # Loop through minModuleSize values
  for (i in seq_along(minModuleSizes)) {
    minSize <- minModuleSizes[i]
    clusters <- dynamicTreeCut::cutreeDynamicTree(hc_fit,
                                                  deepSplit = FALSE,
                                                  minModuleSize = minSize)

    # Calculate metrics
    nClusters[i] <- length(unique(clusters))
  }

  # Calculate first-order differences of nClusters and
  # find the index of the largest negative difference
  max_diff_index <- which.min(diff(nClusters))

  # Correct the index to match the original minModuleSizes
  corrected_index <- max_diff_index + 1

  # Return the optimal minModuleSize
  return(minModuleSizes[corrected_index])
}
