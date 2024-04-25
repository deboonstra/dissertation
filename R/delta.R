# The function of this script file is to define a function to calculate a
# measure of disparity between two matrices, as a way to select the proper
# correlation structure.
delta <- function(a, b, penalty = 0, log = FALSE) {
  # Checking parameter values
  if (!("matrix" %in% class(a))) {
    stop("a must be a matrix.")
  }
  if (!("matrix" %in% class(b))) {
    stop("b must be a matrix.")
  }
  if (all(dim(a) != dim(b)) && (dim(a)[1] != dim(a)[2])) {
    stop("a and b must be square matrices of the same dimension.")
  }
  if (!is.numeric(penalty)) {
    stop("penalty must be a numeric.")
  }
  if (!is.logical(log)) {
    stop("log must be a logical.")
  }

  # Calculating measure ####

  ## Getting diagonal elments ####
  a_diag <- diag(a)
  b_diag <- diag(b)

  ## Calculating the two forms of matrix products ####
  ainv_b <- b_diag / a_diag
  binv_a <- a_diag / b_diag

  ### log-transforming ratios
  if (log) {
    ainv_b <- log(ainv_b)
    binv_a <- log(binv_a)
  }

  ## Calculating measure ####
  measure <- sum(pmax(ainv_b, binv_a)) + penalty

  # Returning value ####
  return(measure)
}