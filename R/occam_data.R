# The function of this script file is to generate a responses based solely on
# N, n, and p for any distribution.

occam_data <- function(N, n, p, family = stats::gaussian(), beta = NULL) {
  # checking parameter types
  if (!(N %% 1 == 0) || !(n %% 1 == 0) || !(p %% 1 == 0)) {
    stop("N, n, and p must be of integer type.")
  }

  if (class(family) != "family") {
    stop("family must be of family type.")
  }
  valid_dist <- c("gaussian", "poisson", "binomial", "Gamma")
  if (!(family$family %in% valid_dist)) {
    stop("Distribution specified cannot be simulated.")
  }

  if (!is.null(beta)) {
    if (!is.numeric(beta) || length(beta) != p) {
      stop("beta must be numeric and of the same length as p.")
    }
  }

  # Generating covariate matrix ####

  ## Initializing covariate matrix
  x <- matrix(
    data = NA, nrow = N * n, ncol = p,
    dimnames = list(
      NULL,
      c("intercept", paste0("X", seq_len(p - 1)))
    )
  )

  ## Intercept term ####
  x[, 1] <- rep(x = 1, times = nrow(x))

  ## Remaining covariates ####
  for (j in 2:p) {
    x[, j] <- runif(n = nrow(x), min = -1, max = 1)
  }
  # Generating mean structure parameters ####
  if (is.null(beta)) {
    beta <- rep(x = 1, times = p)
  } else {
    beta <- beta
  }

  # Generating response vector ####

  ## Calculating eta ####
  eta <- c(family$linkinv(x %*% beta))

  ## Calculating response vector ####
  if (family$family == "gaussian") {
    y <- rnorm(n = length(eta), mean = eta, sd = 1)
  } else if (family$family == "poisson") {
    y <- rpois(n = length(eta), lambda = eta)
  } else if (family$family == "binomial") {
    y <- rbinom(n = length(eta), size = 1, prob = eta)
  } else if (family$family == "Gamma") {
    y <- rgamma(n = length(eta), shape = eta, rate = 1)
  }

  # Creating id vector ####
  id <- rep(seq_len(N), each = n)

  # Returning simulated data ###
  out <- list(y = y, X = x, id = id, family = family)
  return(out)
}
