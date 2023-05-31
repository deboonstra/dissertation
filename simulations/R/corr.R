# The function of this script file is to be a collection of R functions used
# to create correlation structures

# ar1 ####
ar1 <- function(n, rho) {
    x <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
    rho^x
}

# cs (compound symmetric/exchangeable) ####
cs <- function(n, rho) {
    mat <- diag(nrow = n, ncol = n)
    mat[upper.tri(mat) | lower.tri(mat)] <- rho
    mat
}

# indp ####
indp <- function(n) {
    diag(nrow = n, ncol = n)
}