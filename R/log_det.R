log_det <- function(object, modified = FALSE, env = parent.frame()) {
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

  # Calculating KL divergence

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

  ## Calculating KL divergence
  log_det <- log(base::det(ai) / base::det(vr))

  # Output
  return(log_det)
}