set.seed(42)
n <- 1000
p <- 1500
beta_cardinality <- 10


data <- ctlars::generate_ccg_data(n,
                                  p,
                                  mean_xr = 0,
                                  sd_xr = 1,
                                  mean_xi = 0,
                                  sd_xi = 1,
                                  beta_cardinality = beta_cardinality,
                                  noise_power = 1,
                                  set_snr = TRUE,
                                  snr_is_linear = TRUE,
                                  snr_val_linear = 2,
                                  snr_val_db = 10
                                  )
cat("True support: ", data$support, "\n")


ctrex_res <- CTrex::ctrex(data$X,
      data$y,
      tFDR = 0.1,
      K = 20,
      max_num_dummies = 10,
      max_T_stop = TRUE,
      method = "ctrex",
      parallel_process = FALSE,
      parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
      seed = NULL,
      eps = .Machine$double.eps,
      verbose = TRUE)

cat("True support: ", sort(data$support), "\n")
cat("Estimated support: ", sort(which(ctrex_res$selected_var > 0)), "\n")

cat("FDP: ", ctlars::FDP(sort(which(ctrex_res$selected_var > 0)), data$support), "\n")
cat("TPP: ", ctlars::TPP(sort(which(ctrex_res$selected_var > 0)), data$support), "\n")
