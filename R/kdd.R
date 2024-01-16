# The function of this script file is to create a new version of
# Kullback-Leilber divergence equation for parameter estimates

kdd <- function(
  object1, object2, symmetric = FALSE, mean_included = FALSE, modified = FALSE,
  env = parent.frame()
) {
  # checking object type
  if (!("geeglm" %in% class(object1))) {
    stop("object must be a geepack::geeglm object.")
  }
  if (!("geeglm" %in% class(object2))) {
    stop("object must be a geepack::geeglm object.")
  }
  if (!is.logical(symmetric)) {
    stop(
      paste(
        "symmetric must be a logical denoting if symmetric divergence is to",
        "be calculated."
      )
    )
  }
  if (!is.logical(mean_included)) {
    stop(
      paste(
        "mean_included must be a logical denoting if the mean difference",
        "should be included in the calculation."
      )
    )
  }
  if (length(object1$geese$beta) != length(object2$geese$beta)) {
    stop("Dimensionality of object1 and object2 do NOT match.")
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

  # Calculating KDD

  ## Obtaining covariance matrices
  vr1 <- object1$geese$vbeta
  if (!modified) {
    object2$call$corstr <- "independence"
    model_indep <- eval(object2$call, envir = env)
    mb2 <- model_indep$geese$vbeta.naiv
  } else {
    mb2 <- object2$geese$vbeta.naiv
  }

  ## Obtaining dimension of mean parameter estimates
  k <- length(object1$geese$beta)

  ## Obtaining mean parameter estimates
  mu1 <- matrix(data = object1$geese$beta, nrow = k, ncol = 1, byrow = TRUE)
  mu2 <- matrix(data = object2$geese$beta, nrow = k, ncol = 1, byrow = TRUE)

  ## Obtaining the components of KDD

  ### Covariance components
  log_det <- log(det(mb2) / det(vr1))
  tr <- sum(diag(solve(mb2) %*% vr1))

  ### Mean difference component
  if (mean_included) {
    mean_diff <- crossprod(x = mu2 - mu1, y = solve(mb2)) %*% (mu2 - mu1)
  } else {
    mean_diff <- 0
  }

  ## Obtaining KDD
  kdd <- 0.5 * (log_det + tr + mean_diff - k)

  ## Obtaining KSD if symmetric is true
  if (symmetric) {
    tr_flip <- sum(diag(solve(vr1) %*% mb2))
    if (mean_included) {
      mean_diff_flip <- crossprod(x = mu1 - mu2, y = solve(vr1)) %*% (mu1 - mu2)
    } else {
      mean_diff_flip <- 0
    }
    kl <- kdd + 0.5 * ((1 / log_det) + tr_flip + mean_diff_flip - k)
  } else {
    kl <- kdd
  }

  # Output
  return(c(kl))
}