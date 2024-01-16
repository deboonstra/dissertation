# The function of this script file is to create a function called bqicu that
# produces the BIC version of QICu. Additionally, this script file cleans up
# how information criteria is pulled from the QIC function given by geepack.

# quasilik ####
quasilik <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputting quasi-likelihood
  unname(geepack::QIC(object)[3])
}

# nparms ####
nparms <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputting number of parameters
  unname(geepack::QIC(object)[5])
}

# gof ####
gof <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputing goodness-of-fit
  -2 * quasilik(object)
}

# qic ####
qic <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputting qic
  unname(geepack::QIC(object)[1])
}

# qicu ####
qicu <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputting qicu
  unname(geepack::QIC(object)[2])
}

# cic ####
cic <- function(object1, object2, modified = FALSE, env = parent.frame()) {
  # checking object type ####
  if (!("geeglm" %in% class(object1))) {
    stop("object1 must be a geepack::geeglm object.")
  }
  if (!missing(object2)) {
    if (!("geeglm" %in% class(object2))) {
      stop("object2 must be a geepack::geeglm object.")
    }
    if (length(object1$geese$beta) != length(object2$geese$beta)) {
      stop("Dimensionality of object1 and object2 do NOT match.")
    }
  }
  if (!is.logical(modified)) {
    stop("modified must be a logical value.")
  }
  if (class(env) != "environment") {
    stop("env must be an environment object.")
  }

  # calculating cic ####
  vr <- object1$geese$vbeta
  if (missing(object2)) {
    if (!modified) {
      cic <- unname(geepack::QIC(object1)[4])
    } else {
      omega <- solve(object1$geese$vbeta.naiv)
      cic <- sum(diag(omega %*% vr))
    }
  } else {
    if (!modified) {
      object2$call$corstr <- "independence"
      model_indep <- eval(object2$call, envir = env)
      omega <- solve(model_indep$geese$vbeta.naiv)
    } else {
      omega <- solve(object2$geese$vbeta.naiv)
    }
    cic <- sum(diag(omega %*% vr))
  }

  # returning output ####
  return(cic)
}

# qicc ####
qicc <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }
  # outputting ####
  unname(geepack::QIC(object)[6])
}

# bqicu ####
bqicu <- function(object) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }

  # calculating bqicu
  n <- length(unique(object$id))
  gof(object) + (log(n) * nparms(object))
}