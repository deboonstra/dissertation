# The function of this script file is to run a simulation to determine the
# effectiveness of using the GLS approach with GEE to better determine the
# GOF term in QIC. See the 2023_04_20 and 2023_05_18 notes for more detail.

# Loading libraries and functions ####
R <- list.files(path = "./simulations/R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5

# Simulation ####
set.seed(1997)
res <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
    ## Simulating data ####
    dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
    ## Transform y ####
    yy <- transform_y(dat)
    ## Produce best subsets
    best_sub <- expand_cols(dat$X)
    # Fitting independence model
    qic <- rep(NA, length(best_sub))
    for (i in seq_along(best_sub)) {
        XX <- dat$X[, best_sub[[i]]]
        f1 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "independence")
        qic[i] <- unname(geepack::QIC(f1)[1])
    }
    # Determining the best_sub index ####
    # that resulted in the minimum QIC value
    w <- which.min(qic)
    # Updating res object ####
    res[[j]] <- list(
        qic = qic, min = w, qic_min = qic[w], vars = best_sub[[w]]
    )
    setTxtProgressBar(pb, j)
}
close(pb)

# Save simulations ####
saveRDS(res, file = "./simulations/outputs/norm_gls_sim.rds")