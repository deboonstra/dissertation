# The function of this script fileis to create a function named sim_data2 that
# simulates a normal response vector using a method different from sim_data,
# which is based on the Wang et al (2012) paper. The response vector will be
# generated from design matrix that includes a intercept, one binary explanatory
# variable, and p - 2 variables generated from the uniform distribution.
sim_norm <- function(
  N, n, beta, rho = 0.5, corstr = "exchangeable", zcor = NULL
) {

  # checking parameter types
  if (!(N %% 1 == 0) || !(n %% 1 == 0)) {
    stop("N and n must be of integer type.")
  }

  if (length(beta) == 0) {
    stop("beta must be at least a length of 1.")
  }

  if ((rho < -1) || (rho > 1)) {
    stop("rho must be between -1 and 1.")
  }

  if (!(corstr %in% c("exchangeable", "ar1", "independence", "userdefined"))) {
    stop(
      paste0(
        "corstr must be one of the following: exchangeable, ar1,",
        "independence, userdefined."
      )
    )
  }

  if (corstr == "userdefined") {
    if (is.null(zcor)) {
      stop("zcor must be specified if corstr = userdefined")
    }

    if (!is.null(zcor) && !is.matrix(zcor)) {
      stop("zcor must be a square matrix")
    }

    if ((dim(zcor)[1] != n) || (dim(zcor)[2] != n)) {
      stop("zcor must be a square matrix of size n x n.")
    }
  }

  p <- length(beta)
  # initializing covariate matrix
  tmpts <- rep(1:n, times = N)
  X <- matrix(NA, nrow = N * n, ncol = p)

  # generating continuous covariates
  # assigning values of covariates
  X[, 1] <- rep(x = 1, times = N * n) # intercept
  X[, 2] <- rbinom(n = N * n, size = 1, prob = 0.5) # binary variable

  ## continuous variables
  even <- which(seq_len(p) %% 2 == 0)
  for (j in 3:p) {
    if (j %in% even) {
      X[, j] <- stats::runif(length(tmpts), tmpts, tmpts + 1)
    } else {
      X[, j] <- -1 * stats::runif(length(tmpts), tmpts, tmpts + 1)
    }
  }

  # generating within-subject effect
  if (corstr == "exchangeable") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = cs(n = n, rho = rho)
    )
  } else if (corstr == "ar1") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = ar1(n = n, rho = rho)
    )
  } else if (corstr == "independence") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = indp(n = n)
    )
  } else if (corstr == "userdefined") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = zcor
    )
  }
  epsilon <- c(t(epsilon))

  # creating id vector
  id <- rep(seq_len(N), each = n)

  # generating response vector
  y <- c(X %*% beta + epsilon)

  # returning structued list of class "normal"
  colnames(X) <- c("Intercept", paste0("X", seq_len(p - 1)))
  names(beta) <- c("Intercept", paste0("X", seq_len(p - 1)))
  return(
    structure(list(y = y, X = X, id = id, beta = beta), class = "sim.data")
  )
}