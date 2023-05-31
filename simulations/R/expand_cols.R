# The function of this script file is to create a function that produces a
# list containing all the best subset combinations based on a covariate matrix.
# This function does not return the subsetted covariate matrix but the column
# indecies.

expand_cols <- function(X) {
    p <- ncol(X)
    unlist(
        sapply(seq_len(p), function(x) combn(p, x, simplify = FALSE)),
        recursive = FALSE
    )
}