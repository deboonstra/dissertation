# The function of this script file is to be a collection of script files used
# to simulate data. Currently, this is only written to generate a normal
# response based on the wang_2012 paper (PGEEs). For other reponse types,
# a `dist` parameter and if-then(s) may easily be added.

sim_data <- function(N, n, beta, rho) {
    p <- length(beta)
    # initializing covariate matrix
    X <- matrix(NA, nrow = N * n, ncol = p)
    # generating continuous covariates
    mat <- mvtnorm::rmvnorm(
        n = N * (p - 1),
        mean = rep(0, n),
        sigma = ar1(n = n, rho = 0.5)
    )
    mat <- c(t(mat))
    mat <- matrix(mat, nrow = N * n, ncol = (p - 1))
    # assigning values of covariates
    X[, 1] <- rbinom(N * n, 1, 0.5)
    X[, 2:p] <- mat
    # generating within-subject effect
    epsilon <- mvtnorm::rmvnorm(
        n = N, mean = rep(0, n), sigma = cs(n = n, rho = rho)
    )
    epsilon <- c(t(epsilon))
    # creating id vector
    id <- rep(seq_len(N), each = n)
    # generating response vector
    y <- X %*% beta + epsilon
    # returning structued list of class "normal"
    colnames(X) <- paste0("X", seq_len(p))
    names(beta) <- paste0("X", seq_len(p))
    return(structure(list(y = y, X = X, id = id, beta = beta), class = "normal"))
}