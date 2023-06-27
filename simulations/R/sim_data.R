# The function of this script file is to be a collection of script files used
# to simulate data. Currently, this is only written to generate a normal
# response based on the wang_2012 paper (PGEEs). For other reponse types,
# a `dist` parameter and if-then(s) may easily be added.

sim_data <- function(
  N, n, beta, family = stats::gaussian(), rho = 0.5,
  corstr = "exchangeable", zcor = NULL) {

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

  if ("family" != class(family)) {
    stop("family must be a family object produce by stats::family.")
  }

  if (!(corstr %in% c("exchangeable", "ar1", "independence", "userdefined"))) {
    stop("corstr must be one of the following: exchangeable, ar1, independence, userdefined.")
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
  X <- matrix(NA, nrow = N * n, ncol = p)
  # generating continuous covariates
  mat <- mvtnorm::rmvnorm(
    n = N * (p - 1),
    mean = rep(0, n),
    sigma = ar1(n = n, rho = 0.5)
  )
  mat <- c(t(mat))
  mat <- matrix(mat, nrow = N * n, ncol = (p - 1))
  # assigning values of covariates
  X[, 1] <- rbinom(N * n, 1, 0.5)
  X[, 2:p] <- mat
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
  # fam <- glm(1 ~ 1, family = family)$family #nolint
  # eta <- fam$linkfun(c(X %*% beta + epsilon)) #nolint
  # y <- rdist(N * n, eta) #nolint
  y <- c(X %*% beta + epsilon)
  # returning structued list of class "normal"
  colnames(X) <- paste0("X", seq_len(p))
  names(beta) <- paste0("X", seq_len(p))
  return(structure(list(y = y, X = X, id = id, beta = beta), class = "normal"))
}