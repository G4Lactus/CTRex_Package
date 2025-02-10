library(cmvnorm)

# Non-circulant symmetric case
# with nnd covariance matrix
# with non-zero pseudo/relation matrix
# --------------------------------------------------------------------------
set.seed(123)

p <- 4
n <- 1e4

# Create a Hermitian positive semi-definite cov mat
A <- matrix(complex(real = stats::rnorm(p**2),
                     imaginary = stats::rnorm(p**2)),
             nrow = p, ncol = p)
Sigma <- Conj(t(A)) %*% A
diag(Sigma) <- 1
cat("\nCovariance Matrix:\n")
print(Sigma)

# Generate samples from complex circularly symmetric distribution
samples <- cmvnorm::rcmvnorm(n, mean = rep(0, times = p), sigma = Sigma)
head(samples)
plot(Re(samples), Im(samples), col = "blue")
grid()
corrplot::corrplot(cor(Re(samples), Im(samples)))


# Check empirical covariance
scaled_samples <- scale(samples, center = TRUE, scale = FALSE)
emp_cov <- Conj(t(scaled_samples)) %*% scaled_samples  / (n - 1)
cat("\nEmpirical Covariance Matrix:\n")
print(emp_cov)

emp_pcov <- t(scaled_samples) %*% scaled_samples / (n - 1)
cat("\nEmpirical Pseudo-Covariance Matrix:\n")
print(emp_pcov)


# Compute variances (diagonal elements)
vars <- diag(emp_cov)

# Compute the correlation matrix
corr_mat <- emp_cov / sqrt(outer(vars, vars, "*"))

# Visualization
# correlation matrix based on covariance matrix
corrplot::corrplot(Mod(corr_mat), method = "color", title = "Magnitude cormat")
corrplot::corrplot(Re(corr_mat), method = "color", title = "Real part cormat")
corrplot::corrplot(Im(corr_mat), method = "color", title = "Imag part cormat")
corrplot::corrplot(cor(Re(corr_mat), Im(corr_mat)),
                   method = "color",
                   title = "Correlation Re and Im")


# pseudo-covariance matrix
corrplot::corrplot(Mod(emp_pcov), method = "color", title = "Magnitude Pseudo-Covariance Matrix")
corrplot::corrplot(Re(emp_pcov), method = "color", title = "Magnitude Pseudo-Covariance Matrix")
corrplot::corrplot(Im(emp_pcov), method = "color", title = "Magnitude Pseudo-Covariance Matrix")
corrplot::corrplot(cor(Re(emp_pcov), Im(emp_pcov)),
                   method = "color",
                   title = "Correlation Re and Im")

