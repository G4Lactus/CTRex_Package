n <- 100
p <- 300

X <- matrix(data = complex(real = stats::rnorm(n * p),
                           imaginary = stats::rnorm(n * p)),
            nrow = n, ncol = p)

X[, 25:50] <- matrix(data = )

phases <- Arg(X)

# Compute angular distances for clustering
angular_distance <- function(theta1, theta2) {
  diff <- abs(theta1 - theta2)
  return(min(diff, 2 * pi - diff))
}

# Example usage for clustering
dist_matrix <- matrix(0, n, n)
for (i in 1:n) {
  for (j in i:n) {
    theta1 <- phases[i, ]
    theta2 <- phases[j, ]
    dist_matrix[i, j] <- mean(sapply(1:p, function(k) angular_distance(theta1[k], theta2[k])))
    dist_matrix[j, i] <- dist_matrix[i, j]
  }
}

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
plot(hc, main = "Clustering Based on Phase Features")
