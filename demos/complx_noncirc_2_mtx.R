generate_complex_vector <- function(n,
                                    real_mean = 0,
                                    real_cov = NULL,
                                    imag_mean = 0,
                                    imag_cov = NULL) {

  # Set default covariance matrices if not provided
  if (is.null(real_cov)) real_cov <- diag(length(real_mean))
  if (is.null(imag_cov)) imag_cov <- diag(length(imag_mean))

  # Generate real and imaginary parts
  real_part <- mvtnorm::rmvnorm(n, mean = real_mean, sigma = real_cov)
  imag_part <- mvtnorm::rmvnorm(n, mean = imag_mean, sigma = imag_cov)

  # Combine into complex vector
  complex_vector <- complex(real = real_part, imaginary = imag_part)

  return(complex_vector)
}


library(mvtnorm)

# Improper Case
n <- 100
real_mean <- c(0, 1)
imag_mean <- c(-1, 0)
real_cov <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
imag_cov <- matrix(c(1, -0.3, -0.3, 1), nrow = 2)

improp_samples <- generate_complex_vector(n,
                                   real_mean,
                                   real_cov,
                                   imag_mean,
                                   imag_cov)

plot(Mod(improp_samples), type="l")
grid()
plot(Arg(improp_samples), type="l")
grid()


# Proper Case
prop_samples <- generate_complex_vector(n,
                                          real_mean,
                                          real_cov,
                                          real_mean,
                                          real_cov)

plot(Mod(prop_samples), type="l")
grid()
plot(Arg(prop_samples), type="l")
grid()




generate_proper_complex_normal <- function(n_samples, mean, cov) {
  # Ensure mean is a complex vector
  if (!is.complex(mean)) stop("Mean must be a complex vector")

  # Dimension of the complex vector
  p <- length(mean)

  # Augmented covariance matrix
  augmented_cov <- rbind(
    cbind(cov, matrix(0, nrow = p, ncol = p)),
    cbind(matrix(0, nrow = p, ncol = p), cov)
  )

  # Generate augmented real samples
  augmented_samples <- MASS::mvrnorm(n_samples, mu = rep(0, 2 * p), Sigma = augmented_cov)

  # Split into real and imaginary parts
  real_part <- augmented_samples[, 1:p] + Re(mean)
  imag_part <- augmented_samples[, (p + 1):(2 * p)] + Im(mean)

  # Combine into complex samples
  complex_samples <- real_part + 1i * imag_part

  return(complex_samples)
}

# Load required libraries
library(MASS)   # For mvrnorm
library(ggplot2) # For visualization

# Define parameters
n_samples <- 1000
mean <- c(1 + 2i, -1 + 1i)    # Complex mean vector
cov <- matrix(c(1, 0.5,       # Covariance matrix for real and imaginary parts
                0.5, 2),
              nrow = 2)

# Generate samples
set.seed(42) # For reproducibility
Z <- generate_proper_complex_normal(n_samples, mean, cov)

# Extract real and imaginary parts
real_part <- Re(Z)
imag_part <- Im(Z)

# Create a data frame for visualization
data <- data.frame(
  Real = as.vector(real_part),
  Imaginary = as.vector(imag_part)
)

# Scatter plot of real vs imaginary parts
ggplot(data, aes(x = Real, y = Imaginary)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Complex Normal Distribution Samples",
       x = "Real Part",
       y = "Imaginary Part")

# Density plot of real and imaginary parts
ggplot(data) +
  geom_density(aes(x = Real), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = Imaginary), fill = "orange", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density of Real and Imaginary Parts",
       x = "Value",
       y = "Density")



generate_improper_complex_normal <- function(n_samples, mean, cov, pseudo_cov) {
  # Ensure mean is a complex vector
  if (!is.complex(mean)) stop("Mean must be a complex vector")

  # Dimension of the complex vector
  p <- length(mean)

  # Construct augmented covariance matrix
  augmented_cov <- rbind(
    cbind(Re(cov + pseudo_cov), Im(pseudo_cov - 1i * cov)),
    cbind(Im(pseudo_cov + 1i * cov), Re(cov - pseudo_cov))
  )

  # Generate augmented real samples
  augmented_samples <- MASS::mvrnorm(n_samples, mu = rep(0, 2 * p), Sigma = augmented_cov)

  # Split into real and imaginary parts
  real_part <- augmented_samples[, 1:p] + Re(mean)
  imag_part <- augmented_samples[, (p + 1):(2 * p)] + Im(mean)

  # Combine into complex samples
  complex_samples <- real_part + 1i * imag_part

  return(complex_samples)
}


# Parameters for improper distribution
n_samples <- 1000
mean <- c(1 + 2i, -1 + 1i)    # Complex mean vector
cov <- matrix(c(1.0, 0.5+2i,       # Covariance matrix for real and imaginary parts
                0.5-2i, 2.0),
              nrow = 2)

pseudo_cov <- matrix(c(0.3, -0.1 + 0.4i,   # Pseudo-covariance matrix
                       -0.1 + 0.4i, -0.2),
                     nrow = 2)

# Generate improper samples
set.seed(42)
Z_improper <- generate_improper_complex_normal(n_samples, mean, cov, pseudo_cov)

# Extract real and imaginary parts
real_part_improper <- Re(Z_improper)
imag_part_improper <- Im(Z_improper)

# Create a data frame for visualization
data_improper <- data.frame(
  Real = as.vector(real_part_improper),
  Imaginary = as.vector(imag_part_improper)
)

# Scatter plot of real vs imaginary parts (improper case)
library(ggplot2)
ggplot(data_improper, aes(x = Real, y = Imaginary)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Improper Complex Normal Distribution Samples",
       x = "Real Part",
       y = "Imaginary Part")
