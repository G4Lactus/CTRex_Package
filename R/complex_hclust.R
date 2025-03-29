#' @title
#' Perform Hierarchical Clustering on Complex-Valued Data
#'
#' @description
#' Internal function that conducts hierarchical clustering on complex-valued
#' predictors using correlation-based distances. Supports both fixed-height and
#' dynamic tree cutting methods for cluster identification.
#' Critical component for Group-Variable Selection (GVS) in CT-Rex methodology.
#'
#' @param X Complex-valued predictor matrix (n x p).
#' @param hc_method Hierarchical clustering method: \itemize{
#'   \item "single": Single linkage
#'   \item "complete": Complete linkage
#'   \item "ward.D2": Ward's method
#' }
#' @param hc_tree_cut_type Dendrogram segmentation method: \itemize{
#'   \item "fixed": Cut at specified height
#'   \item "dynamic": Dynamic tree cutting
#' }
#' @param hc_tree_cut_height Height threshold for fixed cutting (0-1).
#' Default: 0.5.
#' @param hc_dyn_tree_minModuleSize Minimum cluster size for dynamic cutting.
#' Default: 10.
#'
#' @return List containing: \itemize{
#'   \item max_clusters: Maximum number of identified clusters
#'   \item clusters: Data frame with variable cluster assignments
#'   \item cluster_sizes: Vector of cluster member counts
#' }
#'
#' @importFrom stats hclust cutree aggregate
#' @importFrom dynamicTreeCut cutreeDynamicTree
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' X <- matrix(complex(real=rnorm(100), imaginary=rnorm(100)), ncol=20)
#'
#' # Fixed height cutting
#' clust_fixed <- complex_hclust(X, "ward.D2", "fixed", hc_tree_cut_height=0.6)
#'
#' # Dynamic cutting (requires dynamicTreeCut package)
#' if(requireNamespace("dynamicTreeCut", quietly = TRUE)) {
#'   clust_dynamic <- complex_hclust(X, "complete", "dynamic")
#' }
#'
#' }
#
complex_hclust <- function(X,
                           hc_method,
                           hc_tree_cut_type,
                           hc_tree_cut_height = 0.5,
                           hc_dyn_tree_minModuleSize = 10) {

  # TODO: Input data fidelity check

  # Perform hierarchical clustering
  hc_fit <- stats::hclust(complex_corr_dist(X), method = hc_method)

  if (hc_tree_cut_type == "fixed") {
    clusters <- stats::cutree(hc_fit, h = 1 - hc_tree_cut_height)
    max_clusters <- max(clusters)

  } else if (hc_tree_cut_type == "dynamic") {
    if (is.null(hc_dyn_tree_minModuleSize)) {
      # Perform model order determination
      minModuleSizes <- seq(5, 50, by = 5)
      optimal_minModuleSize <- determine_optimal_minModuleSize(hc_fit, minModuleSizes)
      hc_dyn_tree_minModuleSize <- optimal_minModuleSize

    }
    clusters <- dynamicTreeCut::cutreeDynamicTree(
                        hc_fit,
                        deepSplit = FALSE,
                        minModuleSize = hc_dyn_tree_minModuleSize
                      )
    clusters <- clusters + 1
    names(clusters) <- paste0("V", 1:length(clusters))
    max_clusters <- max(clusters)

  } else {
    stop("`hc_tree_cut_type` not supported.")
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



#' @title
#' Determine Optimal Minimum Cluster Size
#'
#' @description
#' Internal helper function that automatically determines the optimal minimum
#' cluster size for dynamic tree cutting by analyzing cluster count stability
#' across different size thresholds.
#'
#' @param hc_fit Hierarchical clustering object from stats::hclust
#' @param minModuleSizes Vector of candidate minimum cluster sizes to evaluate
#'
#' @return Optimal minimum cluster size from input candidates
#'
#' @keywords internal
#' @importFrom dynamicTreeCut cutreeDynamicTree
#
determine_optimal_minModuleSize <- function(hc_fit, minModuleSizes) {
  # Initialize vectors to store metrics
  nClusters <- numeric(length(minModuleSizes))

  # Loop through hc_dyn_tree_minModuleSize values
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

  # Return the optimal hc_dyn_tree_minModuleSize
  return(minModuleSizes[corrected_index])
}
