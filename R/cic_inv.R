cic_inv <- function(object, modified = FALSE, env = parent.frame()) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }

  if (!is.logical(modified)) {
    stop(
      paste(
        "modified must be logical denoting if model-based working",
        "covariance should be used."
      )
    )
  }

  if (class(env) != "environment") {
    stop("env must be an environment object.")
  }

  # Calculating inverted CIC

  ## Obtaining covariance matrices
  if (!modified) {
    object$call$corstr <- "independence"
    model_indep <- eval(object$call, envir = env)
    vr <- object$geese$vbeta
    ai <- model_indep$geese$vbeta.naiv
  } else {
    vr <- object$geese$vbeta
    ai <- object$geese$vbeta.naiv
  }

  ## Inverted CIC
  cic_inv <- sum(diag(solve(vr) %*% ai))

  # Output
  return(cic_inv)
}