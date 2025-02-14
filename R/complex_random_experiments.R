#' @title
#' Run K random experiments
#'
#' @description
#' Run K random experiments and compute the matrix of relative occurrences for all variables
#' and all numbers of included variables before stopping.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param K Number of random experiments.
#' @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
#' @param num_dummies Number of dummies that are appended to the predictor matrix.
#' @param method 'ctrex' for the T-Rex selector. 'ctrex+GVS' for CT-Rex Group-Variable Selector.
#' @param gvs_type 'cEN' for elastic network group selection. 'cIEN' for informed-elastic network group selection.
#' @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included.
#' Otherwise the entire solution path is computed.
# @param lars_state_list If parallel_process = TRUE: List of state variables of the previous T-LARS steps of the K random experiments
# (necessary for warm-starts, i.e., restarting the forward selection process exactly where it was previously terminated).
# If parallel_process = FALSE: List of objects of the class tlars_cpp associated with the K random experiments
# (necessary for warm-starts, i.e., restarting the forward selection process exactly where it was previously terminated).
#' @param ctlars_learner_lst A list of ctlars learners.
#' @param verbose Logical. If TRUE progress in computations is shown.
# @param intercept Logical. If TRUE an intercept is included.
#' @param intercept Default FALSE.
#' @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
#' @param parallel_process Logical. If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing
#' (default: minimum{Number of random experiments K, number of physical cores}).
#' @param seed Seed for random number generator (ignored if parallel_process = FALSE).
#' @param eps Numerical zero.
#'
#' @return List containing the results of the K random experiments.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParWorkers registerDoSEQ `%do%` `%dopar%` foreach
#' @importFrom doRNG `%dorng%`
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' set.seed(123)
# ----------------------------------------------------------------------
complex_random_experiments <- function(X,
                               y,
                               K = 20,
                               T_stop = 1,
                               num_dummies = ncol(X),
                               dummy_type = c("Complex Gaussian", "Complex Circular Uniform")[1],
                               method = c("ctrex", "ctrex+GVS")[1],
                               gvs_type = c("cEN", "cIEN")[1],
                               hc_method = "ward.D2",
                               coherence_max = 0.5,
                               lambda_2_lars = NULL,
                               early_stop = TRUE,
                               ctlars_learner_lst = NULL,
                               verbose = TRUE,
                               intercept = FALSE,
                               standardize = TRUE,
                               parallel_process = FALSE,
                               parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                               seed = NULL,
                               eps = .Machine$double.eps) {

    # TODO: data fidelity check
    # ...

    # Setup complex Lars learner
    n <- length(y)
    p <- ncol(X)
    y_copy <- y

    res <- list()
    for (k in seq.int(K)) {

      if (T_stop == 1) { # ----------------------------------------------------

        if (method == "ctrex") { # ---------------------------------------------

          # TODO: create environment based wrapper
          # -------------------------------------------
          # Create dummies
          if (dummy_type == "Complex Gaussian") {
            # CN(0, 1) dummies
            X_dummy <- matrix(
              data = complex(
                real = stats::rnorm(n*num_dummies, mean = 0, sd = 1),
                imaginary = stats::rnorm(n*num_dummies, mean = 0, sd = 1)
              )/sqrt(2),
              nrow = n, ncol = num_dummies)
          } else if (dummy_type == "Complex Circular Uniform") {
            # exp(1i U(0, 2pi)) dummies
            X_dummy <- matrix(
              data = exp(1i * stats::runif(n*num_dummies, min = 0, max = 2*pi)),
              nrow = n, ncol = num_dummies)
          } else {
            stop("Dummy type not supported.")
          }
          X_dummy <- cbind(X, X_dummy)
          # -------------------------------------------


        } else if (method == "ctrex+GVS") { # ----------------------------------

          y <- y_copy

          # Ridge regression to determine lambda_2 for elastic net
          if (is.null(lambda_2_lars)) {
            lambda_2_lars <- lambda_cv(X,
                                       y,
                                       epsilon = 1e-3,
                                       L = 100,
                                       nfolds = 10)$lambda_1se
          }

          # trex+GVS data modification
          # -----------------------------------
          # Complex Elastic network
          if (gvs_type == "cEN") { # ----------------------------------

            X_dummy <- add_complex_gvs_dummies(
              X = X,
              num_dummies = num_dummies,
              gvs_type = gvs_type,
              hc_method = hc_method,
              coherence_max = coherence_max)$X_dummy

            p_dummy <- ncol(X_dummy)

            X_dummy <- (1 / sqrt(1 + lambda_2_lars)) *
              rbind(X_dummy, diag(rep(sqrt(lambda_2_lars), times = p_dummy)))

            y <- append(y, rep(0, times = p_dummy))


          # Complex Informed Elastic Network
          } else if (gvs_type == "cIEN") { # ----------------------------------

            cmplx_gvs_dummies <- add_complex_gvs_dummies(
              X = X,
              num_dummies = num_dummies,
              gvs_type = gvs_type,
              hc_method = hc_method,
              coherence_max = coherence_max)

            X_dummy <- cmplx_gvs_dummies$X_dummy

            p_dummy <- ncol(X_dummy)

            max_clusters <- cmplx_gvs_dummies$max_clusters

            cluster_sizes <- cmplx_gvs_dummies$cluster_sizes

            IEN_cl_id_vectors <- cmplx_gvs_dummies$IEN_cl_id_vectors

            X_dummy <- sqrt(lambda_2_lars) *
              rbind((1 / sqrt(lambda_2_lars)) * X_dummy,
                              (1 / sqrt(cluster_sizes)) *
                      matrix(rep(IEN_cl_id_vectors, times = p_dummy / p),
                             ncol = ncol(IEN_cl_id_vectors) * p_dummy / p))

            y <- append(y, rep(0, times = max_clusters))

          } else { # ----------------------------------------------------------
            stop("`gvs_type` not supported.")
          }

          # Re-scale X_dummy and y
          X_dummy <- scale_x(X_dummy)
          y <- y - mean(y)

        } else { # ------------------------------------------------------------
          stop("`method` not supported.")
        }


        # Perform random experiments
        # ----------------------------------------------
        # create learner
        ctlars_learner <- ctlars::ctlars$new(
          X_dummy,
          y,
          standardize = TRUE,
          has_intercept = FALSE,
          num_dummies = num_dummies,
          verbose = verbose
        )

        # conduct first T_step
        ctlars_learner$execute_clars_step(
          t_stop = T_stop,
          early_stop = TRUE,
          use_chol = TRUE
        )
        # cat("Aktive set of learner ", k, ": ", ctlars_learner$get_active_set(), "\n")

        # summarize ctlars learner results
        ctlars_learner_lst[[k]] <- ctlars_learner


      } else { # T_stop > 1 ---------------------------------------------------
        ctlars_learner_lst[[k]]$execute_clars_step(
          t_stop = T_stop,
          early_stop = TRUE,
          use_chol = TRUE
        )
        # cat("Aktive set of learner ", k, ": ", ctlars_learner_lst[[k]]$get_active_set(), "\n")
      }

      # Number of included dummies along solution path
      lars_path <- ctlars_learner_lst[[k]]$get_beta_history()

      # dummy detector
      dummy_num_path <- colSums(matrix(
        abs(lars_path[seq.int(from = p+1, to = p+num_dummies), ]) > eps,
        nrow = num_dummies,
        ncol = ncol(lars_path)
      ))

      # Relative occurrences
      phi_T_mat <- matrix(0, nrow = p, ncol = T_stop)
      for (c in seq(T_stop)) {
        if (!any(dummy_num_path == c)) {
          ind_sol_path <- length(dummy_num_path)
          warning(
            paste(
              "For T_stop = ",
              c,
              " LARS is running until k = min(n, p) and stops there before selecting ",
              c,
              " dummies.",
              sep = ""
            )
          )
        } else {
          ind_sol_path <- which(as.numeric(dummy_num_path) == c)[1]
        }
        phi_T_mat[, c] <- (1 / K) * (abs(lars_path[1:p, ind_sol_path]) > eps)
      }

      model_name <- paste0("mod_", k)
      res[[model_name]]$phi_T_mat <- phi_T_mat
      res[[model_name]]$rand_exp_last_betas_mat <- lars_path[1:p, ncol(lars_path)]
  }

  phi_T_mat <- Reduce("+", lapply(res, function(resX) { return(resX$phi_T_mat) }))
  rand_exp_last_betas_mat <- unname(Reduce(rbind, lapply(res, function(resX) {
    return( resX$rand_exp_last_betas_mat ) })))
  Phi <- apply(abs(rand_exp_last_betas_mat) > eps, 2, sum) / K

  # List of results
  # -----------------------
  rand_exp_res <- list(
    phi_T_mat = phi_T_mat,
    rand_exp_last_betas_mat = rand_exp_last_betas_mat,
    Phi = Phi,
    ctlars_learner_lst = ctlars_learner_lst,
    K = K,
    T_stop = T_stop,
    num_dummies = num_dummies,
    method = method,
    seed = seed,
    eps = eps
  )

  return(rand_exp_res)
}





data_fidelity_check <- function() {
  # # Error control
  # method <- match.arg(method, c("trex", "trex+GVS"))
  #
  # if (!is.matrix(X)) {
  #   stop("'X' must be a matrix.")
  # }
  #
  # if (!is.complex(X)) {
  #   stop("'X' only allows numerical values.")
  # }
  #
  # if (anyNA(X)) {
  #   stop("'X' contains NAs. Please remove or impute them before proceeding.")
  # }
  #
  # if (!is.vector(drop(y))) {
  #   stop("'y' must be a vector.")
  # }
  #
  # if (!is.complex(y)) {
  #   stop("'y' only allows numerical values.")
  # }
  #
  # if (anyNA(y)) {
  #   stop("'y' contains NAs. Please remove or impute them before proceeding.")
  # }
  #
  # if (nrow(X) != length(drop(y))) {
  #   stop("Number of rows in X does not match length of y.")
  # }
  #
  # if (length(K) != 1 ||
  #     K < 2 ||
  #     K %% 1 != 0) {
  #   stop("The number of random experiments 'K' must be an integer larger or equal to 2.")
  # }
  #
  # if (method == "trex") {
  #   if (length(num_dummies) != 1 ||
  #       num_dummies %% 1 != 0 ||
  #       num_dummies < 1) {
  #     stop("'num_dummies' must be an integer larger or equal to 1.")
  #   }
  # }
  #
  # # Number of variables in X
  # p <- ncol(X)
  #
  # # Continue error control
  # if (method == "trex+GVS") {
  #   if (length(num_dummies) != 1 ||
  #       num_dummies %% p != 0 ||
  #       num_dummies < 1) {
  #     stop(
  #       "`num_dummies` must be a positive integer multiple of the total number of original predictors in X."
  #     )
  #   }
  # }
  #
  # if (length(T_stop) != 1 || !(T_stop %in% seq(1, num_dummies))) {
  #   stop(
  #     paste0(
  #       "Value of 'T_stop' not valid. 'T_stop' must be an integer from 1 to ",
  #       num_dummies,
  #       "."
  #     )
  #   )
  # }
  #
  # if (method == "trex+GVS") {
  #   if (length(corr_max) != 1 ||
  #       corr_max < 0 ||
  #       corr_max > 1) {
  #     stop("'corr_max' must have a value between zero and one.")
  #   }
  #
  #   if (!is.null(lambda_2_lars)) {
  #     if (length(lambda_2_lars) != 1 || lambda_2_lars < eps) {
  #       stop("'lambda_2_lars' must be a number larger than zero.")
  #     }
  #   }
  # }
  #
  # if (parallel_process &&
  #     (
  #       length(parallel_max_cores) != 1 ||
  #       parallel_max_cores %% 1 != 0 ||
  #       parallel_max_cores < 2
  #     )) {
  #   stop(
  #     "For parallel processing at least two workers have to be registered:
  #        'parallel_max_cores' must be an integer larger or equal to 2."
  #   )
  # }
  #
  # if (parallel_process &&
  #     parallel_max_cores > min(K, max(1, parallel::detectCores(logical = FALSE)))) {
  #   parallel_max_cores_modified <-
  #     min(K, max(1, parallel::detectCores(logical = FALSE)))
  #   message(
  #     paste0(
  #       "For computing ",
  #       K,
  #       " random experiments, it is not useful/possible to register ",
  #       parallel_max_cores,
  #       " workers. Setting parallel_max_cores = ",
  #       min(K, max(
  #         1, parallel::detectCores(logical = FALSE)
  #       )),
  #       " (# physical cores) ...\n"
  #     )
  #   )
  #   parallel_max_cores <- min(K, max(1, parallel::detectCores(logical = FALSE)))
  # }
  #
  # if (parallel_process && T_stop == 1 && num_dummies <= p) {
  #   message(
  #     "Computing random experiments in parallel...
  #           Note that this is only advantageous if you have at least a few thousand predictors and/or data points in 'X'.
  #           Otherwise, the overhead will slow down the computations in parallel. Thus, for small data sizes it is better
  #           to set parallel_process = FALSE.
  #
  #     Be careful!"
  #   )
  # }
  #
  # if (!(missing(lars_state_list) || is.null(lars_state_list))) {
  #   if (length(lars_state_list) != K) {
  #     stop("Length of 'lars_state_list' must be equal to number of random experiments 'K'.")
  #   }
  # }
}
