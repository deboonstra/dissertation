# The function of this script file is to create a function that calculates
# Kullback-Leibler divergence using gee objects. This calculation is based on
# the asymptotic distribution of the mean structure parameter estimates being
# multivariate normal.

# object, gee object
# symmetric, logical if symmetric divergence should be calculated
# tol, tolerance for inverting matrix

kd <- function(object, symmetric = FALSE, tol = .Machine$double.eps) {
  # checking parameter types ####
  if (!("gee" %in% class(object))) {
    stop("object must be a gee::gee object.")
  }
  if (!is.logical(symmetric)) {
    stop(
      paste(
        "symmetric must be a logical denoting if symmetric divergence is to",
        "be calculated."
      )
    )
  }
  if (!(class(tol) %in% c("numeric", "integer"))) {
    stop("tol must be a numeric or an integer value")
  }

  # Setup functions
  invert <- if ("MASS" %in% loadedNamespaces()) {
    MASS::ginv
  } else {
    solve
  }

  # Calculating Kullback divergence ####

  ## Calculating components of divergence ####

  ###  Number mean parameter estimates ####
  k <- length(object$coefficients)

  ### Covariance matrices ####
  mb <- object$naive.variance
  vr <- object$robust.variance

  ### Calculating trace terms ####
  tr <- sum(diag(crossprod(x = invert(mb, tol = tol), y = vr)))
  if (symmetric) {
    tr_flip <- sum(diag(crossprod(x = invert(vr, tol = tol), y = mb)))
  }

  ### log-determinant terms ####
  ld <- log(det(mb) / det(vr))

  ## Calculating divergence ####
  ## The log-determinant components cancel out for symmetric divergence
  if (symmetric) {
    kl <- 0.5 * (tr + tr_flip) - k
  } else {
    kl <- 0.5 * (ld + tr - k)
  }

  # Output
  return(kl)
}