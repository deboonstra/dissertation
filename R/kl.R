kl <- function(
  object, type = "kdd", modified = FALSE, env = parent.frame()
) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }

  type <- tolower(type)
  if (!(type %in% c("kdd", "ksd"))) {
    stop("type must be a character value denoting either KDD or KSD.")
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
  log_diff_kdd <- log(base::det(ai) / base::det(vr))
  kdd <- 0.5 * (log_diff_kdd + sum(diag(solve(ai) %*% vr)) - ncol(vr))

  if (type == "kdd") {
    ### Calculating KDD
    kl <- kdd
  } else {
    ### Calculating KSD
    log_diff_ksd <- log(base::det(vr) / base::det(ai))
    kl <- kdd + 0.5 * (log_diff_ksd + sum(diag(solve(vr) %*% ai)) - ncol(ai))
  }

  # Output
  return(kl)
}