# The function of this script file is to create a function called bqicu that
# produces the BIC version of QICu. Additionally, this script file cleans up
# how information criteria is pulled from the QIC function given by geepack.

# quasilik ####
quasilik <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputting quasi-likelihood
    unname(geepack::QIC(object)[3])
}

# nparms ####
nparms <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputting number of parameters
    unname(geepack::QIC(object)[5])
}

# gof ####
gof <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputing goodness-of-fit
    -2 * quasilik(object)
}

# qic ####
qic <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputting qic
    unname(geepack::QIC(object)[1])
}

# qicu ####
qicu <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputting qicu
    unname(geepack::QIC(object)[2])
}

# cic ####
cic <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputing cic
    unname(geepack::QIC(object)[4])
}

# qicc ####
qicc <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }
    # outputting ####
    unname(geepack::QIC(object)[6])
}

# bqicu ####
bqicu <- function(object) {
    # checking object type
    if ("geeglm" %notin% class(object)) {
        stop("bqicu requires a geeglm object as input.")
    }

    # calculating bqicu
    n <- length(unique(object$id))
    gof(object) + (log(n) * nparms(object))
}