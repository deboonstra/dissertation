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
res_indp <- list()
res_cs <- list()
res_ar1 <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
    ## Simulating data ####
    dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
    ## Transform y ####
    yy <- transform_y(dat)
    ## Produce best subsets
    best_sub <- expand_cols(dat$X)
    # Fitting independence model
    qic1 <- rep(NA, length(best_sub))
    qic2 <- rep(NA, length(best_sub))
    qic3 <- rep(NA, length(best_sub))
    for (i in seq_along(best_sub)) {
        XX <- dat$X[, best_sub[[i]]]
        f1 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "independence")
        f2 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "exchangeable")
        f3 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "ar1")
        qic1[i] <- unname(geepack::QIC(f1)[1])
        qic2[i] <- unname(geepack::QIC(f2)[1])
        qic3[i] <- unname(geepack::QIC(f3)[1])
    }
    # Determining the best_sub index ####
    # that resulted in the minimum QIC value
    w1 <- which.min(qic1)
    w2 <- which.min(qic2)
    w3 <- which.min(qic3)
    # Updating res object ####
    res_indp[[j]] <- list(
        qic = qic1, min = w1, qic_min = qic1[w1], vars = best_sub[[w1]]
    )
    res_cs[[j]] <- list(
        qic = qic2, min = w2, qic_min = qic2[w2], vars = best_sub[[w2]]
    )
    res_ar1[[j]] <- list(
        qic = qic3, min = w3, qic_min = qic3[w3], vars = best_sub[[w3]]
    )
    setTxtProgressBar(pb, j)
}
close(pb)

# Save simulations ####
save(res_indp, res_cs, res_ar1, file = "./simulations/outputs/norm_gls_sim.RData")