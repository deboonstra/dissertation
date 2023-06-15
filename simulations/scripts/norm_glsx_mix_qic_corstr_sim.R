# The function of this script file is to run a simulation to determine the
# effectiveness of using a mixed GLS approach with GEEs. The GOF term will be
# a by-product of the transformed y and X; however, the penalty term will be
# based on the original data.

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
    y <- dat$y
    X <- dat$X
    id <- dat$id
    ## Transform data ####
    transform_data <- transform_y(dat)
    yt <- transform_data$y
    Xt <- transform_data$X
    ## Produce best subsets ####
    best_sub <- expand_cols(X)
    ## Fitting models ####
    ic <- gf <- pen <- corstr <- rep(NA, length(best_sub))
    wc <- c("independence", "exchangeable", "ar1")
    for (i in seq_along(best_sub)) {
        XX <- X[, best_sub[[i]]]
        XXt <- Xt[, best_sub[[i]]]
        f1 <- geepack::geeglm(y ~ XX, id = id, corstr = wc[1])
        f2 <- geepack::geeglm(y ~ XX, id = id, corstr = wc[2])
        f3 <- geepack::geeglm(y ~ XX, id = id, corstr = wc[3])
        ft <- geepack::geeglm(yt ~ XXt, id = id, corstr = "independence")
        ### Selecting best performing corstr ####
        ww <- which.min(c(cic(f1), cic(f2), cic(f3)))
        gf[i] <- gof(ft)
        pen[i] <- 2 * c(cic(f1), cic(f2), cic(f3))[ww]
        corstr[i] <- wc[ww]
        ic[i] <- gf[i] + pen[i]
    }
    ## Determining the best_sub index ####
    ## that resulted in the minimum QIC
    w <- which.min(ic)
    ## Updating res object ####
    res[[j]] <- list(
        qic = ic, gof = gf, penalty = pen, corstr = corstr,
        min = w, qic_min = ic[w], corstr_min = corstr[w], vars = best_sub[[w]]
    )
    setTxtProgressBar(pb, j)
    if (j == nsims) {
        close(pb)
    }
}

# Save simulation ####
if (!dir.exists("./simulations/outputs/norm_glsx_mix_qic_corstr_sim/")) {
    dir.create("./simulations/outputs/norm_glsx_mix_qic_corstr_sim/")
}
saveRDS(
    res,
    file = "./simulations/outputs/norm_glsx_mix_qic_corstr_sim/norm_glsx_mix_qic_corstr_sim.rds"
)
