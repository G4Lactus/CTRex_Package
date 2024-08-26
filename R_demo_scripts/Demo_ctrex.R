n <- 1000
p <- 15000
beta_cardinality <- 10
set.seed(42)

data <- ctlars::generate_ccg_data(n,
                                  p,
                                  mean_xr = 0,
                                  sd_xr = 1,
                                  mean_xi = 0,
                                  sd_xi = 1,
                                  beta_cardinality = beta_cardinality,
                                  noise_power = 1,
                                  set_snr = FALSE,
                                  snr_is_linear = TRUE,
                                  snr_val_linear = 2,
                                  snr_val_db = 10
                                  )
cat("True support: ", data$support, "\n")


ctrex_res <- ctrex(data$X,
      data$y,
      tFDR = 0.1,
      K = 20,
      max_num_dummies = 10,
      max_T_stop = TRUE,
      method = "trex",
      parallel_process = FALSE,
      parallel_max_cores = min(K, max(1, parallel::detectCores(logical = FALSE))),
      seed = NULL,
      eps = .Machine$double.eps,
      verbose = TRUE)

cat("True support: ", sort(data$support), "\n")
cat("Estimated support: ", sort(which(ctrex_res$selected_var > 0)))

cat("FDP: ", FDP(sort(which(ctrex_res$selected_var > 0)), data$support), "\n")
cat("TPP: ", TPP(sort(which(ctrex_res$selected_var > 0)), data$support), "\n")
