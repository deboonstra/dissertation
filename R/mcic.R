mcic <- function(object, tol = .Machine$double.eps) {
  # Checking function parameter types ####
  if (!("gee" %in% class(object))) {
    stop("object must be an object produced by gee::gee")
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

  # Obtaining covariance matrices ####

  ## Inverse of model-based covariance under independence ####
  omega <- invert(object$naive.variance, tol = tol)

  ## Robust covariance given working correlation structure ####
  vr <- object$robust.variance

  # Calculating mCIC ####
  trace <- sum(diag(x = crossprod(x = omega, y = vr)))
  names(trace) <- "mCIC"

  # Returning mCIC
  return(trace)
}