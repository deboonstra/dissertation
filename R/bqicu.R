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
cic <- function(object, modified = FALSE) {
  # checking object type
  if (!("geeglm" %in% class(object))) {
    stop("object must be a geepack::geeglm object.")
  }

  if (!is.logical(modified)) {
    stop("modified must be a logical value.")
  }
  
  # outputing cic
  if (!modified) {
    cic <- unname(geepack::QIC(object)[4])
  } else {
    vr <- object$geese$vbeta
    omega <- solve(object$geese$vbeta.naiv)
    cic <- sum(diag(omega %*% vr))
  }
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