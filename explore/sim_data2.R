# The function of this script is to create a function that simulates data for
# any exponential family distribution supported by the stats package.

sim_data2 <- function(
  N, n, beta, family = stats::gaussian(), rho = 0.5,
  corstr = "exchangeable", zvar = 1, zcor = NULL) {

  # checking parameter values
  if (!(N %% 1 == 0) || !(n %% 1 == 0)) {
    stop("N and n must be of integer type.")
  }

  if (length(beta) == 0) {
    stop("beta must be at least a length of 1.")
  }

  if ((rho < -1) || (rho > 1)) {
    stop("rho must be between -1 and 1.")
  }

  if (zvar < 0 || length(zvar) > 1) {
    stop("zvar must be a scalar greater than 0.")
  }

  if ("family" != class(family)) {
    stop("family must be a family object produce by stats::family.")
  }

  if (!(corstr %in% c("exchangeable", "ar1", "independence", "userdefined"))) {
    stop(
      paste0(
        "corstr must be one of the following: ",
        "exchangeable, ar1, independence, userdefined."
      )
    )
  }

  if (corstr == "userdefined") {
    if (is.null(zcor)) {
      stop("zcor must be specified if corstr = userdefined")
    }

    if (!is.null(zcor) && !is.matrix(zcor)) {
      stop("zcor must be a square p.s.d matrix")
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
  X[, 1] <- rbinom(N * n, 1, 0.5)
  for (j in 2:p) {
    X[, j] <- log(stats::runif(length(tmpts), tmpts, tmpts + 1))
  }
  # generating within-subject effect
  if (corstr == "exchangeable") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = zvar * cs(n = n, rho = rho)
    )
  } else if (corstr == "ar1") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = zvar * ar1(n = n, rho = rho)
    )
  } else if (corstr == "independence") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = zvar * indp(n = n)
    )
  } else if (corstr == "userdefined") {
    epsilon <- mvtnorm::rmvnorm(
      n = N, mean = rep(0, n), sigma = zvar * zcor
    )
  }
  epsilon <- c(t(epsilon))
  # creating id vector
  id <- rep(seq_len(N), each = n)
  # generating response vector
  fam <- glm(1 ~ 1, family = family)$family
  eta <- fam$linkinv(c(X %*% beta + epsilon))
  if (fam$family == "gaussian") {
    y <- eta
  } else if (fam$family == "poisson") {
    y <- stats::rpois(N * n, eta)
  } else if (fam$family == "binomial") {
    y <- stats::rbinom(N * n, 1, eta)
  } else if (fam$family == "Gamma") {
    y <- stats::rgamma(N * n, eta)
  }
  # returning structued list of class "normal"
  colnames(X) <- paste0("X", seq_len(p))
  names(beta) <- paste0("X", seq_len(p))
  return(structure(list(y = y, X = X, id = id, beta = beta), class = "sim.data")) #nolint
}