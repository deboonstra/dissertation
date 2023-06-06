# The function of this script file is to run a simulation study where the
# performance of AIC and QIC will be compared. This comparison can only occur if the working correlation structure is independence.

# Loading libraries and functions ####
R <- list.files(path = "./simulations/R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0
corstr <- "independence"

# Simulation ####
set.seed(1997)
res_qic <- list()
res_aic <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
    ## Simulating data ####
    dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr)
    ## Produce best subsets ####
    best_sub <- expand_cols(dat$X)
    # Fitting models ####
    q <- gof1 <- penalty1 <- rep(NA, length(best_sub))
    aic <- gof2 <- penalty2 <- rep(NA, length(best_sub))
    for (i in seq_along(best_sub)) {
        XX <- dat$X[, best_sub[[i]]]
        f1 <- geepack::geeglm(dat$y ~ XX, id = dat$id, corstr = corstr)
        f2 <- stats::glm(dat$y ~ XX)
        q[i] <- qic(f1)
        gof1[i] <- gof(f1)
        penalty1[i] <- cic(f1)
        aic[i] <- stats::AIC(f2)
        gof2[i] <- -2 * logLik(f2)
        penalty2[i] <- aic[i] - gof2[i]
    }
    # Determining the best_sub index ####
    # that resulted in the minimum AIC/QIC values
    w1 <- which.min(q)
    w2 <- which.min(aic)
    # Updating res objects ####
    res_qic[[j]] <- list(
        qic = q, gof = gof1, penalty = penalty1,
        min = w1, qic_min = q[w1], vars = best_sub[[w1]]
    )
    res_aic[[j]] <- list(
        aic = aic, gof = gof2, penalty = penalty2,
        min = w2, qic_min = aic[w2], vars = best_sub[[w2]]
    )
    setTxtProgressBar(pb, j)
    if (j == nsims) {
        close(pb)
    }
}

# Save simulations ####
save(
    res_qic, res_aic,
    file = "./simulations/outputs/norm_aic_sim/norm_aic_sim.RData"
)