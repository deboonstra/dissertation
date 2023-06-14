# The function of this script file is to run a simulation study, where the 
# performance of QIC, QICu, and BQICu will be compared through the traditional
# application of GEEs and information criterions (i.e., no transformations).
# We are wanting to investigate these criterions performance in selecting the
# correct mean and correlation structure. norm_bqicu_sim only focused on
# selecting the proper mean structure.

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
res_qic <- list()
res_qicu <- list()
res_bqicu <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
    ## Simulatin data ####
    dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
    y <- dat$y
    X <- dat$X
    id <- dat$id
    ## Producing best subsets
    best_sub <- expand_cols(dat$X)
    ## Fitting models ####
    ## with different correlation structures
    Q <- gfQ <- penQ <- corQ <- rep(NA, length(best_sub))
    u <- gfu <- penu <- coru <- rep(NA, length(best_sub))
    b <- gfb <- penb <- corb <- rep(NA, length(best_sub))
    corstr <- c("independence", "exchangeable", "ar1")
    for (i in seq_along(best_sub)) {
        XX <- X[, best_sub[[i]]]
        f1 <- geepack::geeglm(y ~ XX, id = id, corstr = corstr[1])
        f2 <- geepack::geeglm(y ~ XX, id = id, corstr = corstr[2])
        f3 <- geepack::geeglm(y ~ XX, id = id, corstr = corstr[3])
        ### Selecting best performing corstr ####
        ### based on each IC
        f <- list(f1, f2, f3)
        QQ <- c(qic(f1), qic(f2), qic(f3))
        uu <- c(qicu(f1), qicu(f2), qicu(f3))
        bb <- c(bqicu(f1), bqicu(f2), bqicu(f3))
        wQ <- which.min(QQ)
        wu <- which.min(uu)
        wb <- which.min(bb)
        Q[i] <- QQ[wQ]
        gfQ[i] <- gof(f[[wQ]])
        penQ[i] <- 2 * cic(f[[wQ]])
        corQ[i] <- corstr[wQ]
        u[i] <- uu[wu]
        gfu[i] <- gof(f[[wu]])
        penu[i] <- u[i] - gfu[i]
        coru[i] <- corstr[wu]
        b[i] <- bb[wu]
        gfb[i] <- gof(f[[wb]])
        penb[i] <- b[i] - gfb[i]
        corb[i] <- corstr[wb]
    }
    ## Determining the best_sub index ####
    ## that resulted in the minimum QIC,
    ## QICu, and BQICu values
    wwQ <- which.min(Q)
    wwu <- which.min(u)
    wwb <- which.min(b)
    ## Updating res object(s) ####
    res_qic[[j]] <- list(
        qic = Q, gof = gfQ, penalty = penQ, corstr = corQ,
        min = wwQ, qic_min = Q[wwQ], corstr_min = corQ[wwQ],
        vars = best_sub[[wwQ]]
    )
    res_qicu[[j]] <- list(
        qicu = u, gof = gfu, penalty = penu, corstr = coru,
        min = wwu, qicu_min = u[wwu], corstr_min = coru[wwu],
        vars = best_sub[[wwu]]
    )
    res_bqicu[[j]] <- list(
        bqicu = b, gof = gfb, penalty = penb, corstr = corb,
        min = wwb, bqicu_min = b[wwb], corstr_min = corb[wwb],
        vars = best_sub[[wwb]]
    )
    ## Updating progress ####
    setTxtProgressBar(pb, j)
    if (j == nsims) {
        close(pb)
    }
}

# Save simulations ####
save(
    res_qic, res_qicu, res_bqicu,
    file = "./simulations/outputs/norm_bqicu_corstr_sim/norm_bqicu_corstr_sim.RData"
)