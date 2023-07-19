# initializing covariate matrix
tmpts <- rep(1:n, times = N)
direction <- sample(x = c(1, -1), size = p - 1, replace = TRUE)
X <- matrix(NA, nrow = N * n, ncol = p)
# generating continuous covariates
# assigning values of covariates
X[, 1] <- rbinom(N * n, 1, 0.5)
for (j in 2:p) {
  if (direction[j - 1] == 1) {
    X[, j] <- log(stats::runif(length(tmpts), tmpts, tmpts + 1))
  } else {
    X[, j] <- stats::runif(length(tmpts), 1 / (tmpts + 1), 1 / tmpts)
  }
}