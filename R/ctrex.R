#' @title
#' Run the CT-Rex selector
#'
#' @description
#' The CT-Rex selector performs fast variable selection in high-dimensional settings while
#' controlling the false discovery rate (FDR) at a user-defined target level.
#'
#' @param X Real valued predictor matrix.
#' @param y Response vector.
#' @param tFDR Target FDR level (between 0 and 1, i.e., 0% and 100%).
#' @param K Number of random experiments.
#' @param max_num_dummies Integer factor determining the maximum number of dummies as a multiple of the number of original variables p
#' (i.e., num_dummies = max_num_dummies * p).
#' @param max_T_stop If TRUE the maximum number of dummies that can be included before stopping is set to ceiling(n / 2),
#' where n is the number of data points/observations.
#' @param method 'trex' for the T-Rex selector.
#' @param parallel_process Logical. If TRUE random experiments are executed in parallel.
#' @param parallel_max_cores Maximum number of cores to be used for parallel processing
#' (default: minimum{Number of random experiments K, number of physical cores}).
#' @param seed Seed for random number generator (ignored if parallel_process = FALSE).
#' @param eps Numerical zero.
#' @param verbose Logical. If TRUE progress in computations is shown.
#'
#' @return A list containing the estimated support vector and additional information, including the number of used dummies and the number of included dummies before stopping.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParWorkers registerDoSEQ
#'
#' @export
#'
#' @examples
#' set.seed(1)
# ----------------------------------------------------------------------
ctrex <- function(X,
                 y,
                 tFDR = 0.2,
                 K = 20,
                 max_num_dummies = 10,
                 max_T_stop = TRUE,
                 method = "trex",
                 parallel_process = FALSE,
                 parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
                 seed = NULL,
                 eps = .Machine$double.eps,
                 verbose = TRUE) {

  # ---------------------------------------------------
  # TODO: perform `data_fidelity_check` function
  # ---------------------------------------------------
  # Scale X and center y
  X <- scale_x(X)
  y <- y - mean(y)


  # --------------------------------------------------
  # TODO: move to function `determine_voting_grid`
  # Voting level grid
  V <- seq(0.5, 1 - eps, by = 1 / K)
  V_len <- length(V)

  # Initialize L-loop
  # --------------------------------------------------
  n <- nrow(X)
  p <- ncol(X)

  LL <- 1
  T_stop <- 1
  FDP_hat <- rep(NA, times = V_len)

  # 75% optimization point for determining number of dummies required
  opt_point <- which(abs(V - 0.75) < eps)
  if (length(opt_point) == 0) {
    # If 75% optimization point choose closest optimization point lower than 75%
    opt_point <- length(V[V < 0.75])
  }
  # --------------------------------------------------
  while ((LL <= max_num_dummies &&
          FDP_hat[opt_point] > tFDR) ||
          sum(!is.na(FDP_hat)) == 0) {
    ctlars_learner_lst <- lapply(seq.int(K), function(k) {return(k)})
    num_dummies <- LL * p
    LL <- LL + 1

    # K Random experiments
    suppressWarnings(
      rand_exp <- complex_random_experiments(
        X = X,
        y = y,
        K = K,
        T_stop = T_stop,
        num_dummies = num_dummies,
        method = method,
        early_stop = TRUE,
        ctlars_learner_lst = ctlars_learner_lst,
        verbose = FALSE,
        standardize = TRUE,
        parallel_process = parallel_process,
        parallel_max_cores = parallel_max_cores,
        seed = seed,
        eps = eps
      )
    )

    # Phi_prime
    # ----------------------------
    Phi_prime <- Phi_prime_fun(
      p = p,
      T_stop = T_stop,
      num_dummies = num_dummies,
      phi_T_mat = rand_exp$phi_T_mat,
      Phi = rand_exp$Phi,
      eps = eps
    )

    # FDP_hat
    # ----------------------
    FDP_hat <- fdp_hat(
      V = V,
      Phi = rand_exp$Phi,
      Phi_prime = Phi_prime
    )

    # Print number of included dummies by the extended calibration algorithm.
    if (verbose) {
      cat(paste("\n Appended dummies:", num_dummies, "\n"))
    }
  }

  #browser()
  # Initialize T-loop
  # ------------------------------------------------------
  FDP_hat_mat <- matrix(FDP_hat, nrow = 1)
  Phi_mat <- matrix(rand_exp$Phi, nrow = 1)

  if (max_T_stop) {
    max_T <- min(num_dummies, ceiling(n / 2))
  } else {
    max_T <- num_dummies
  }

  # Reset seed
  if (!is.null(seed)) {
    seed <- seed + 12345
  }

  while (FDP_hat[V_len] <= tFDR && T_stop < max_T) {
    T_stop <- T_stop + 1

    # K Random experiments
    suppressWarnings(
      rand_exp <- complex_random_experiments(
        X = X,
        y = y,
        K = K,
        T_stop = T_stop,
        num_dummies = num_dummies,
        method = method,
        early_stop = TRUE,
        ctlars_learner_lst = rand_exp$ctlars_learner_lst,
        verbose = FALSE,
        intercept = FALSE,
        standardize = TRUE,
        parallel_process = parallel_process,
        parallel_max_cores = parallel_max_cores,
        seed = seed,
        eps = eps
      )
    )


    Phi_mat <- rbind(Phi_mat, rand_exp$Phi)


    # Phi_prime
    Phi_prime <- Phi_prime_fun(
      p = p,
      T_stop = T_stop,
      num_dummies = num_dummies,
      phi_T_mat = rand_exp$phi_T_mat,
      Phi = rand_exp$Phi,
      eps = eps
    )

    # FDP_hat
    FDP_hat <- fdp_hat(
      V = V,
      Phi = rand_exp$Phi,
      Phi_prime = Phi_prime
    )
    FDP_hat_mat <- rbind(FDP_hat_mat, FDP_hat)

    # Print number of included dummies by the extended calibration algorithm.
    if (verbose) {
      cat(paste("\n Included dummies before stopping:", T_stop, "\n"))
    }
  }
  # -----------------------------------------------

  # T-Rex: Select variables
  # -----------------------------------------------
  res_T_dummy <- select_var_fun(
    p = p,
    tFDR = tFDR,
    T_stop = T_stop,
    FDP_hat_mat = FDP_hat_mat,
    Phi_mat = Phi_mat,
    V = V
  )
  # -----------------------------------------------


  # Compile results
  # -----------------------------------------------
  res <- list(
    selected_var = res_T_dummy$selected_var,
    tFDR = tFDR,
    T_stop = T_stop,
    num_dummies = num_dummies,
    V = V,
    v_thresh = res_T_dummy$v_thresh,
    FDP_hat_mat = FDP_hat_mat,
    Phi_mat = Phi_mat,
    R_mat = res_T_dummy$R_mat,
    phi_T_mat = rand_exp$phi_T_mat,
    Phi_prime = Phi_prime,
    method = method
  )
  # -----------------------------------------------

  return(res)
}


# Auxiliar functions
# -----------------------------------------------
data_fidelity_check <- function() {
  # Error control
  method <- match.arg(method, c("trex", "trex+GVS"))

  type <- match.arg(type, c("lar", "lasso"))

  if (!is.matrix(X)) {
    stop("'X' must be a matrix.")
  }

  if (!is.complex(X)) {
    stop("'X' only allows only complex numerical values.")
  }

  if (anyNA(X)) {
    stop("'X' contains NAs. Please remove or impute them before proceeding.")
  }

  if (!is.vector(drop(y))) {
    stop("'y' must be a vector.")
  }

  if (!is.complex(y)) {
    stop("'y' only allows numerical values.")
  }

  if (anyNA(y)) {
    stop("'y' contains NAs. Please remove or impute them before proceeding.")
  }

  if (nrow(X) != length(drop(y))) {
    stop("Number of rows in X does not match length of y.")
  }

  if (length(tFDR) != 1 ||
      tFDR < 0 ||
      tFDR > 1) {
    stop("'tFDR' must be a number between 0 and 1 (including 0 and 1).")
  }

  if (length(K) != 1 ||
      K < 2 ||
      K %% 1 != 0) {
    stop("The number of random experiments 'K' must be an integer larger or equal to 2.")
  }

  if (length(max_num_dummies) != 1 ||
      max_num_dummies < 1 ||
      max_num_dummies %% 1 != 0) {
    stop("'max_num_dummies' must be an integer larger or equal to 1.")
  }

  if (parallel_process &&
      (length(parallel_max_cores) != 1 ||
       parallel_max_cores %% 1 != 0 ||
       parallel_max_cores < 2)) {
    stop(
      "For parallel processing at least two workers have to be registered:
         'parallel_max_cores' must be an integer larger or equal to 2."
    )
  }

  if (parallel_process &&
      parallel_max_cores > min(K, max(
        1, parallel::detectCores(logical = FALSE)
      ))) {
    parallel_max_cores_modified <-
      min(K, max(1, parallel::detectCores(logical = FALSE)))
    message(
      paste0(
        "For computing ",
        K,
        " random experiments, it is not useful/possible to register ",
        parallel_max_cores,
        " workers. Setting parallel_max_cores = ",
        min(K, max(
          1, parallel::detectCores(logical = FALSE)
        )),
        " (# physical cores) ...\n"
      )
    )
    parallel_max_cores <-
      min(K, max(1, parallel::detectCores(logical = FALSE)))
  }


}

determine_voting_grid <- function() {

}
