#' @title
#' Generate Complex Dummy Variables for Group Variable Selection
#'
#' @description
#' Internal function that creates complex-valued dummy variables reflecting the
#' correlation structure of original predictors through hierarchical clustering.
#' Core component of CT-Rex's Group-Variable Selection (GVS) methodology.
#'
#' @param X Complex-valued predictor matrix (n x p)
#' @param num_dummies Number of dummy variables to generate. Typically a
#' multiple of p.
#' @param gvs_type Group variable selection type: \itemize{
#'  \item "cEN": Complex Elastic Net.
#'  \item "cIEN": Complex Informed Elastic Net.
#' }
#' @param hc_method Hierarchical clustering method \itemize{
#'  \item "single": Single linkage
#'  \item "complete": Complete linkage
#'  \item "ward.D2": Ward's method
#' } Default: "single".
#' @param hc_tree_cut_type Dendrogram cutting method: \itemize{
#'  \item "fixed": Static height cut
#'  \item "dynamic": Data-driven adaptive clustering
#' } Default: "fixed".
#' @param hc_tree_cut_height Height threshold for tree cutting (0-1).
#' Default 0.5
#' @param hc_dyn_tree_minModuleSize Minimum cluster size for dynamic cutting.
#' Default 10.
#'
#' @return List containing: \itemize{
#'   \item X_dummy: Augmented matrix [X|D] (n x (p + num_dummies))
#'   \item max_clusters: Maximum number of identified clusters
#'   \item cluster_sizes: Vector of cluster member counts
#'   \item IEN_cl_id_vectors: (Only for cIEN) Binary cluster identity matrix
#' }
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Generate complex data
#' X <- matrix(complex(real=rnorm(100), imaginary=rnorm(100)), ncol=20)
#'
#' # Create EN-style dummies with dynamic clustering
#' en_dummies <- add_complex_gvs_dummies(
#'   X,
#'   num_dummies = 50,
#'   gvs_type = "cEN",
#'   hc_method = "ward.D2",
#'   hc_tree_cut_type = "dynamic"
#' )
#'
#' # Create IEN-style dummies with fixed cutting
#' ien_dummies <- add_complex_gvs_dummies(
#'   X,
#'   num_dummies = 30,
#'   gvs_type = "cIEN",
#'   hc_tree_cut_height = 0.7
#' )
#'
#' }
#'
add_complex_gvs_dummies <- function(X,
                                    num_dummies,
                                    gvs_type,
                                    hc_method,
                                    hc_tree_cut_type,
                                    hc_tree_cut_height = 0.5,
                                    hc_dyn_tree_minModuleSize = 10) {

  # TODO: write function for input data fidelity check (see below - open)

  n <- nrow(X)
  p <- ncol(X)
  colnames(X) <- paste0("V", seq(p))

  # Perform clustering
  cluster_res <- complex_hclust(X,
                                hc_method = hc_method,
                                hc_tree_cut_type = hc_tree_cut_type,
                                hc_tree_cut_height = hc_tree_cut_height,
                                hc_dyn_tree_minModuleSize = hc_dyn_tree_minModuleSize
                               )

  # Generate dummy predictors and append them to original X
  X_dummy <- matrix(NA, nrow = n, ncol = p + num_dummies)
  X_dummy[, seq.int(1, p)] <- X

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
    # Binary cluster identity vectors for T-Rex+GVS+IEN
    # ----------------------------------------------------
    IEN_cl_ind <- lapply(cluster_res$clusters$`clusters$Var`,
                         FUN = function(x) as.numeric(sub(".", "", x)))
    IEN_cl_id_vectors <- t(sapply(seq_along(IEN_cl_ind), FUN = function(x) {
      cl_id_vec <- rep(FALSE, times = ncol(X))
      cl_id_vec[IEN_cl_ind[[x]]] <- TRUE
      cl_id_vec}))

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


# TODO: write function for data fidelity check

