# The function of this script file is to define a function to calculate a
# measure of disparity between two matrices, as a way to select the proper
# correlation structure.
delta <- function(a, b, penalty = 0) {
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

  # Calculating measure ####

  ## Getting diagonal elments ####
  a_diag <- diag(a)
  b_diag <- diag(b)

  ## Calculating the two forms of matrix products ####
  ainv_b <- b_diag / a_diag
  binv_a <- a_diag / b_diag

  ## Calculating measure ####
  measure <- sum(pmax(ainv_b, binv_a)) + penalty

  # Returning value ####
  return(measure)
}