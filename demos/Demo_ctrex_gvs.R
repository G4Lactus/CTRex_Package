# Demo: CT-Rex Group Variable Selection with cEN and cIEN
# ===========================================================
n <- 50
p <- 200
p1 <- 20
grp_size <- 5
num_grps <- p1 / grp_size

X <- matrix(data = complex(real = stats::rnorm(n * p),
                           imaginary = stats::rnorm(n * p)),
            nrow = n,
            ncol = p)

generate_cmplx_grp <- function(phase_function, amplitude) {
  phase <- phase_function(1:n)
  random_phases <- stats::runif(n, 0, stats::runif(1, pi/4, pi/2)) # Generate random phases for each element
  cmplx_vals <- amplitude * exp(1i * (phase + random_phases)) # Apply both deterministic and random phases
  return(cmplx_vals)
}

# Phase functions
linear_phase <- function(x) {
  return(0.1 * x)
}

quadratic_phase <- function(x) {
  return(0.5 * x**2)
}

harmonic_phase <- function(x) {
  return(sin(0.5 * x))
}

constant_phase <- function(x) {
  return(rep(pi / 4, length(x)))
}


for (i in (1:grp_size)) {
  X[, i] <- generate_cmplx_grp(linear_phase, stats::runif(n, 1, 2))
}

for (i in (2*grp_size + 1):(3*grp_size)) {
  X[, i] <- generate_cmplx_grp(quadratic_phase, stats::runif(n, 3, 4))
}

for (i in (4*grp_size + 1):(5*grp_size)) {
  X[, i] <- generate_cmplx_grp(harmonic_phase, stats::runif(n, 0.5, 0.9))
}

for (i in (7*grp_size + 1):(8*grp_size)) {
  X[, i] <- generate_cmplx_grp(constant_phase, stats::runif(n, 2.5, 2.9))
}




beta_cardinality <- 5
beta_support <- sample(c((1:grp_size), (7*grp_size + 1):(8*grp_size)),
                       size = beta_cardinality)
beta <- rep(0, times = p)
beta[beta_support] <- complex(
                        real = stats::runif(beta_cardinality, 0, 2*pi),
                        imaginary = stats::runif(beta_cardinality, 0, 2*pi)
                      )

y <- X[,beta_support] %*% beta[beta_support]
y <- y + complex(real = stats::rnorm(n, mean = 0, sd = 1),
                 imaginary = stats::rnorm(n, mean = 0, sd = 1))/sqrt(2)


plot(X)

plot(X[, c((1:grp_size),
           (2*grp_size + 1):(3*grp_size),
           (4*grp_size + 1):(5*grp_size),
           (7*grp_size + 1):(8*grp_size)
           )])

image(Mod(Conj(t(X)) %*%X))
image(Arg(Conj(t(X)) %*%X))

# =============================================================================
# =============================================================================

# Demo
# -------------------------------------------------

# cTRex+gvs selection: cEN
# ------------------------------
ctrex_gvs_en_tc_static <- CTrex::ctrex(
                              X = X,
                              y = y,
                              tFDR = 0.1,
                              method = "ctrex+GVS",
                              gvs_type = "cEN",
                              hc_method = "single",
                              hc_tree_cut_type = "fixed",
                              hc_tree_cut_height = 0.5
                            )


ctrex_gvs_en_tc_dyn <- CTrex::ctrex(
                          X = X,
                          y = y,
                          tFDR = 0.1,
                          method = "ctrex+GVS",
                          gvs_type = "cEN",
                          hc_method = "ward.D2",
                          hc_tree_cut_type = "dynamic",
                          hc_dyn_tree_minModuleSize = 5
                        )


ctrex_gvs_en_tc_dyn_mos <- CTrex::ctrex(
                          X = X,
                          y = y,
                          tFDR = 0.1,
                          method = "ctrex+GVS",
                          gvs_type = "cEN",
                          hc_method = "ward.D2",
                          hc_tree_cut_type = "dynamic",
                          hc_dyn_tree_minModuleSize = NULL
                        )


# cTRex+gvs selection: cIEN
# ------------------------------
ctrex_gvs_ien_tc_static <- CTrex::ctrex(
                       X = X,
                       y = y,
                       tFDR = 0.1,
                       method = "ctrex+GVS",
                       gvs_type = "cIEN",
                       hc_method = "single",
                       hc_tree_cut_type = "fixed",
                       hc_tree_cut_height = 0.5
                      )


ctrex_gvs_ien_tc_dyn <- CTrex::ctrex(
                      X = X,
                      y = y,
                      tFDR = 0.1,
                      method = "ctrex+GVS",
                      gvs_type = "cIEN",
                      hc_method = "ward.D2",
                      hc_tree_cut_type = "dynamic",
                      hc_dyn_tree_minModuleSize = 5
                      )

ctrex_gvs_ien_tc_dyn_mos <- CTrex::ctrex(
                        X = X,
                        y = y,
                        tFDR = 0.1,
                        method = "ctrex+GVS",
                        gvs_type = "cIEN",
                        hc_method = "ward.D2",
                        hc_tree_cut_type = "dynamic",
                        hc_dyn_tree_minModuleSize = NULL
                      )



# Performance Evaluation
# ------------------------------
tlars::TPP(which(Mod(ctrex_gvs_en_tc_static$selected_var) > 0), beta_support)
tlars::TPP(which(ctrex_gvs_en_tc_dyn$selected_var > 0), beta_support)
tlars::TPP(which(ctrex_gvs_ien_tc_static$selected_var > 0), beta_support)
tlars::TPP(which(ctrex_gvs_ien_tc_dyn$selected_var > 0), beta_support)
tlars::TPP(which(ctrex_gvs_en_tc_dyn_mos$selected_var > 0), beta_support)
tlars::TPP(which(ctrex_gvs_ien_tc_dyn_mos$selected_var > 0), beta_support)


tlars::FDP(which(ctrex_gvs_en_tc_static$selected_var > 0), beta_support)
tlars::FDP(which(ctrex_gvs_en_tc_dyn$selected_var > 0), beta_support)
tlars::FDP(which(ctrex_gvs_ien_tc_static$selected_var > 0), beta_support)
tlars::FDP(which(ctrex_gvs_ien_tc_dyn$selected_var > 0), beta_support)
tlars::FDP(which(ctrex_gvs_en_tc_dyn_mos$selected_var > 0), beta_support)
tlars::FDP(which(ctrex_gvs_ien_tc_dyn_mos$selected_var > 0), beta_support)
