# Demo: CT-Rex Group Variable Selection with cEN and cIEN
# =====================================================

# Function setup
# --------------------------------------------------------------------------
generate_grp_covariance <- function(p_group,
                                    phase_method = "random",
                                    rho_mag = 0.9,
                                    phase_slope = 2*pi/10,
                                    phase_freq = 0.1,
                                    phase_amp = pi/3,
                                    phase_alpha = pi/100)  {

  lags <- seq(0, p_group - 1)
  grp_length <- length(p_group)

  # Phase design
  phases <- switch(phase_method,
                   zero = rep(0, grp_length),
                   constant = stats::runif(1, 0, 2*pi),
                   linear = phase_slope * seq.int(0, (grp_length-1)),
                   harmonic = phase_amp * sin(2*pi*phase_freq * (0:(grp_length-1))),
                   quadratic = phase_alpha * seq.int(0, (grp_length-1))**2,
                   random = stats::runif(grp_length, 0, 2*pi),
                   stop("Invalid `phase_method`. Use 'zero', 'constant',
                        'linear', 'quadratic', 'harmonic', or 'random'.")
  )

  # Create first row of Toeplitz matrix with Hermitian structure
  first_row <- rho_mag ** lags * exp(1i * phases)

  # Build Hermitian symmetric Toeplitz matrix manually
  # NOTE: (grp_cov_mat + Conj(t(grp_cov_mat))) / 2 annihilates imaginary part
  grp_cov_mat <- matrix(0, p_group, p_group)
  for (i in 1:p_group) {
    for (j in i:p_group) {
      k <- abs(i - j)
      val <- first_row[k+1]
      grp_cov_mat[i,j] <- val
      grp_cov_mat[j,i] <- Conj(val)
    }
  }

  # Ensure diagonal unity
  diag(grp_cov_mat) <- 1 + 0i

  return(grp_cov_mat)
}


# --------------------------------------------------------------------------
plot_covariance <- function(Gamma,
                            title,
                            plot_type = c("Re", "Im", "Mod", "Arg", "Cor")[1]) {

  if (plot_type %in% c("Re", "Im", "Mod", "Arg")) {
    plot_idx <- seq.int(nrow(Gamma), 1, by = -1)
    Gamma <- t(Gamma[plot_idx,])
  }

  switch(plot_type,
         Re = {
           # Real part
           image(Re(Gamma), main = title,
                 col = grDevices::colorRampPalette(c("navy", "white", "red"))(100)
           )
         },

         Im = {
           # Imaginary part
           image(Im(Gamma), main = title,
                 col = grDevices::colorRampPalette(c("navy", "white", "red"))(100)
           )
         },

         Mod = {
           # Modulus
           image(Mod(Gamma), main = title,
                 col = grDevices::colorRampPalette(c("navy", "white", "red"))(100)
           )
         },

         Arg = {
           image(Arg(Gamma), main = title,
                 col = grDevices::colorRampPalette(c("navy", "white", "red"))(100)
           )
         },

         Cor = {
           corrplot::corrplot(Mod(Gamma), method = "color")
         },

         {
           stop("`plot_type` is not supported.")
         })
}


# --------------------------------------------------------------------------
create_beta <- function(grp_labels,
                        active_groups,
                        beta_type,
                        theta,
                        grp_sizes,
                        phase_method = "random",
                        phase_slope = 2*pi/10,
                        phase_freq = 0.1,
                        phase_amp = pi/3,
                        phase_alpha = pi/100)  {

  # Data fidelity check
  stopifnot(
    sum(grp_sizes) == length(grp_labels),
    all(active_groups %in% unique(grp_labels)),
    phase_method %in% c("random", "zero",
                        "linear", "harmonic",
                        "constant", "quadratic")
  )

  beta <- rep(0 + 0i, times = sum(grp_sizes))

  for (g in active_groups) {
    grp_idxs <- which(g == grp_labels)
    grp_length <- length(grp_idxs)

    # Phase design patterns
    phases <- switch(phase_method,
               random = stats::runif(grp_length, 0, 2*pi),
               constant = stats::runif(1, 0, 2*pi),
               zero = rep(0, grp_length),
               linear = phase_slope * seq.int(0, (grp_length-1)),
               harmonic = phase_amp * sin(2*pi*phase_freq * (0:(grp_length-1))),
               quadratic = phase_alpha * seq.int(0, (grp_length-1))**2,
               stop("Invalid `phase_method`.")
    )


    switch(beta_type,
           cmplx_uniform = {
             beta[grp_idxs] <- exp(1i * phases)
           },

           cmplx_expo_decay = {
             magnitudes <- theta ** (0:(grp_length-1))
             beta[grp_idxs] <- magnitudes * exp(1i * phases)
           },

           cmplx_double_expo_decay = {
             if(grp_length %% 2 == 0) {
               centers <- c(grp_length/2, grp_length/2 + 1)
               distances <- pmin(abs(seq_len(grp_length) - centers[1]),
                                 abs(seq_len(grp_length) - centers[2]))
             } else {
               center <- (grp_length + 1)/2
               distances <- abs(seq_len(grp_length) - center)
             }
             magnitudes <- theta ** distances
             beta[grp_idxs] <- magnitudes * exp(1i * phases)
           },

           { # default case
             stop("`beta_type` not supported.")
           }
    )
  }

  return(beta)
}


# --------------------------------------------------------------------------
create_toeplitz_cmplx_grouped_data <- function(
    n = 100,
    grp_sizes = c(25, 35, 35, 5),
    grp_means = complex(
      real = stats::rnorm(4),
      imaginary = stats::rnorm(4)
    ),
    num_active_grps = 2,
    rho_mag = 0.9,
    theta = 0.5,
    beta_type = c("cmplx_uniform",
                  "cmplx_expo_decay",
                  "cmplx_double_expo_decay")[3],
    phase_method = "random",
    phase_slope = 2*pi/10,
    phase_freq = 0.1,
    phase_amp = pi/3,
    phase_alpha = pi/100,
    snr_linear = 10,
    visualize = TRUE,
    verbose = TRUE)
{

  # Data Input fidelity checks
  num_grps <- length(grp_sizes)

  # total features
  p_features <- sum(grp_sizes)

  # Group labels
  grp_labels <- rep(seq.int(1, num_grps), times = grp_sizes)
  grp_idxs <- lapply(unique(grp_labels), function(label_x) {
    return(which(label_x == grp_labels))
  })

  # Random active group selection
  active_grp_idxs <- sort(sample(seq(1, num_grps - 1), size = num_active_grps))

  # Generate block diagonal covariance matrix
  block_cov <- matrix(0, nrow = p_features, ncol = p_features)

  # Fill block diagonal
  start_idx <- 1
  for (g in seq(1, num_grps)) {
    grp_size_x <- grp_sizes[g]
    end_idx <- start_idx + grp_size_x - 1
    block_grp_idxs <- seq(start_idx, end_idx)
    block_cov[block_grp_idxs, block_grp_idxs] <-
      generate_grp_covariance(
        p_group = grp_size_x,
        phase_method = phase_method,
        rho_mag = rho_mag,
        phase_slope = phase_slope,
        phase_freq = phase_freq,
        phase_amp = phase_amp,
        phase_alpha = phase_alpha)

    start_idx <- end_idx + 1
  }

  if (all.equal(block_cov, Conj(t(block_cov)))) {
    if (verbose) {
      print("Covariance is Hermitian:")
    }
  } else {
    stop("Covariance is not Hermitian.")
  }
  if (visualize) {
    plot_covariance(block_cov, title = "Real part", plot_type = "Re")
    plot_covariance(block_cov, title = "Imaginary part", plot_type = "Im")
    plot_covariance(block_cov, title = "Modulus part", plot_type = "Mod")
    plot_covariance(block_cov, title = "Arg part", plot_type = "Arg")
    plot_covariance(block_cov, title = "Cormat", plot_type = "Cor")
  }

  # Sample the elliptic distribution according to 1st and 2nd moments
  mu <- rep(0, times = p_features)
  start_idx <- 1
  for (g in seq(1, num_grps)) {
    group_size <- grp_sizes[g]
    mu[seq(start_idx, start_idx + group_size - 1)] <- grp_means[g]
    start_idx <- start_idx + group_size
  }
  X <- cmvnorm::rcmvnorm(n,
                         mean = mu,
                         sigma = block_cov)

  # Define support vector
  beta <- create_beta(grp_labels,
                active_groups = active_grp_idxs,
                beta_type,
                theta,
                grp_sizes,
                phase_method,
                phase_slope,
                phase_freq,
                phase_amp,
                phase_alpha)

  if (visualize) {
    par(mfrow = c(2, 2))
    plot(Re(beta), Im(beta), type = "o", col = "blue")
    grid()
    plot(Mod(beta), type = "o", col = "blue")
    grid()
    plot(Arg(beta), type = "o", col = "blue")
    grid()
    par(mfrow = c(1, 1))
  }

  # Simulate response vector y
  active_inds <- unlist(grp_idxs[active_grp_idxs])
  y_true <- X[, active_inds] %*% beta[active_inds]

  # SNR adjustments
  signal_power <- mean(Mod(y_true)**2)
  noise_power <- signal_power / snr_linear

  # white circularly symmetric noise
  noise_sd <- sqrt(noise_power) / 2
  noise <- stats::rnorm(n, mean = 0, sd = noise_sd)

  # final signal
  y <- y_true + noise

  if (visualize) {
    plot(Mod(y), type = "o", col = "blue")
    grid()
  }

  return(
    list(
      y = y,
      X = X,
      active_grp_indices = active_grp_idxs,
      grp_idxs = grp_idxs
    )
  )
}


# ---------------------------------------------------------------------------
eval_performance <- function(selected_vars, group_idxs, active_groups, verbose = TRUE) {

  # Translation into active groups
  sel_active_grps <- unique(
    sapply(which(selected_vars > 0),
           function(idx_x) {
             return(which(sapply(group_idxs, function(group) idx_x %in% group)))
           }))

  fdp <- 0
  tpp <- 0
  if (!is.null(sel_active_grps)) {
    # FDP performance
    fdp <- tlars::FDP(selected_vars = sel_active_grps,
                      true_actives = active_groups)

    # TPP performance
    tpp <- tlars::TPP(selected_vars = sel_active_grps,
                      true_actives = active_groups)
  }
  if (verbose) {
    print(paste0("FDP: ", fdp, " TPP: ", tpp))
  }

  return(c(fdp, tpp))
}



# =============================================================================
# =============================================================================

# Demo
# -----------------------------
#set.seed(234)

cmplx_grp_data <- create_toeplitz_cmplx_grouped_data(
                      n = 80,
                      beta_type = c("cmplx_uniform",
                                    "cmplx_expo_decay",
                                    "cmplx_double_expo_decay")[1],
                      phase_method = "random",
                      rho_mag = 0.5,
                      phase_slope = 2*pi/10,
                      phase_freq = 0.1,
                      phase_amp = pi/3,
                      phase_alpha = pi/100,
                      visualize = TRUE
                    )

browser()

# cTRex+gvs selection: cEN
cgvs_trex_en <- CTrex::ctrex(
                      X = cmplx_grp_data$X,
                      y = cmplx_grp_data$y,
                      tFDR = 0.1,
                      method = "ctrex+GVS",
                      gvs_type = "cEN")


# cTRex+gvs selection: cIEN
cgvs_trex_ien <- CTrex::ctrex(
                       X = cmplx_grp_data$X,
                       y = cmplx_grp_data$y,
                       tFDR = 0.1,
                       method = "ctrex+GVS",
                       gvs_type = "cIEN")


# Performance Evaluation
eval_performance(cgvs_trex_en$selected_var,
                 group_idxs = cmplx_grp_data$grp_idxs,
                 active_groups = cmplx_grp_data$active_grp_indices)


eval_performance(cgvs_trex_ien$selected_var,
                 group_idxs = cmplx_grp_data$grp_idxs,
                 active_groups = cmplx_grp_data$active_grp_indices)
