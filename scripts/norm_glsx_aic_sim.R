# The function of this script file is to run a simulation to determine the
# effectiveness of using the GLS approach with GEE to transform the response and
# covariate matrix. Then, use AIC and BIC to select the mean structure.
# See 2023_06_08 notes for more information.

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
res_aic <- list()
res_bic <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
    ## Simulating data ####
    dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
    ## Transform y ####
    transform_data <- transform_y(dat)
    yy <- transform_data$y
    xx <- transform_data$X
    ## Produce best subsets
    best_sub <- expand_cols(xx)
    # Fitting independence model ####
    aic <- penalty1 <- rep(NA, length(best_sub))
    bic <- penalty2 <- rep(NA, length(best_sub))
    gof <- rep(NA, length(best_sub))
    for (i in seq_along(best_sub)) {
        XX <- xx[, best_sub[[i]]]
        f <- stats::glm(yy ~ XX)
        aic[i] <- stats::AIC(f)
        bic[i] <- stats::BIC(f)
        gof[i] <- -2 * stats::logLik(f)
        penalty1[i] <- aic[i] - gof[i]
        penalty2[i] <- bic[i] - gof[i]
    }
    # Determining the best_sub index ####
    # that resulted in the minimum AIC and BIC
    w1 <- which.min(aic)
    w2 <- which.min(bic)
    # Updating res object ####
    res_aic[[j]] <- list(
        aic = aic, gof = gof, penalty = penalty1,
        min = w1, aic_min = aic[w1], vars = best_sub[[w1]]
    )
    res_bic[[j]] <- list(
        bic = bic, gof = gof, penalty = penalty2,
        min = w2, bic_min = bic[w2], vars = best_sub[[w2]]
    )
    setTxtProgressBar(pb, j)
    if (j == nsims) {
        close(pb)
    }
}

# Save simulations ####
if (!dir.exists("./simulations/outputs/norm_glsx_aic_sim/")) {
    dir.create("./simulations/outputs/norm_glsx_aic_sim/")
}
save(
    res_aic, res_bic,
    file = "./simulations/outputs/norm_glsx_aic_sim/norm_glsx_aic_sim.RData"
)