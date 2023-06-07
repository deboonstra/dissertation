# The function of this script file is to run a simulation to determine the
# effectiveness of using the GLS approach with GEE, where the GOF term in 
# QIC, QICu and BQICu will be adjusted by producing a new y and X.
# See 2023_05_31 notes for more information.

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
res_indp_qic <- list()
res_cs_qic <- list()
res_ar1_qic <- list()
res_indp_qicu <- list()
res_cs_qicu <- list()
res_ar1_qicu <- list()
res_indp_bqicu <- list()
res_cs_bqicu <- list()
res_ar1_bqicu <- list()
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
    qic1 <- penalty1 <- rep(NA, length(best_sub))
    qic2 <- penalty2 <- rep(NA, length(best_sub))
    qic3 <- penalty3 <- rep(NA, length(best_sub))
    qicu1 <- penaltyu1 <- rep(NA, length(best_sub))
    qicu2 <- penaltyu2 <- rep(NA, length(best_sub))
    qicu3 <- penaltyu3 <- rep(NA, length(best_sub))
    bqicu1 <- penaltyb1 <- rep(NA, length(best_sub))
    bqicu2 <- penaltyb2 <- rep(NA, length(best_sub))
    bqicu3 <- penaltyb3 <- rep(NA, length(best_sub))
    gof1 <- gof2 <- gof3 <-rep(NA, length(best_sub))
    for (i in seq_along(best_sub)) {
        XX <- xx[, best_sub[[i]]]
        f1 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "independence")
        f2 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "exchangeable")
        f3 <- geepack::geeglm(yy ~ XX, id = dat$id, corstr = "ar1")
        qic1[i] <- qic(f1)
        qicu1[i] <- qicu(f1)
        bqicu1[i] <- bqicu(f1)
        gof1[i] <- gof(f1)
        penalty1[i] <- cic(f1)
        penaltyu1[i] <- qicu1[i] - gof1[i]
        penaltyb1[i] <- bqicu1[i] - gof1[i]
        qic2[i] <- qic(f2)
        qicu2[i] <- qicu(f2)
        bqicu2[i] <- bqicu(f2)
        gof2[i] <- gof(f2)
        penalty2[i] <- cic(f2)
        penaltyu2[i] <- qicu2[i] - gof2[i]
        penaltyb2[i] <- bqicu2[i] - gof2[i]
        qic3[i] <- qic(f3)
        qicu3[i] <- qicu(f3)
        bqicu3[i] <- bqicu(f3)
        gof3[i] <- gof(f3)
        penalty3[i] <- cic(f3)
        penaltyu3[i] <- qicu3[i] - gof3[i]
        penaltyb3[i] <- bqicu3[i] - gof3[i]
    }
    # Determining the best_sub index ####
    # that resulted in the minimum QIC,
    # QICu, and BQICu values
    w1 <- which.min(qic1)
    w2 <- which.min(qic2)
    w3 <- which.min(qic3)
    wu1 <- which.min(qicu1)
    wu2 <- which.min(qicu2)
    wu3 <- which.min(qicu3)
    wb1 <- which.min(bqicu1)
    wb2 <- which.min(bqicu2)
    wb3 <- which.min(bqicu3)
    # Updating res object ####
    res_indp_qic[[j]] <- list(
        qic = qic1, gof = gof1, penalty = penalty1,
        min = w1, qic_min = qic1[w1], vars = best_sub[[w1]]
    )
    res_cs_qic[[j]] <- list(
        qic = qic2, gof = gof2, penalty = penalty2,
        min = w2, qic_min = qic2[w2], vars = best_sub[[w2]]
    )
    res_ar1_qic[[j]] <- list(
        qic = qic3, gof = gof3, penalty = penalty3,
        min = w3, qic_min = qic3[w3], vars = best_sub[[w3]]
    )
    res_indp_qicu[[j]] <- list(
        qicu = qicu1, gof = gof1, penalty = penaltyu1,
        min = wu1, qicu_min = qicu1[wu1], vars = best_sub[[wu1]]
    )
    res_cs_qicu[[j]] <- list(
        qicu = qicu2, gof = gof2, penalty = penaltyu2,
        min = wu2, qicu_min = qicu2[wu2], vars = best_sub[[wu2]]
    )
    res_ar1_qicu[[j]] <- list(
        qicu = qicu3, gof = gof3, penalty = penaltyu3,
        min = wu3, qicu_min = qicu3[wu3], vars = best_sub[[wu3]]
    )
    res_indp_bqicu[[j]] <- list(
        bqicu = bqicu1, gof = gof1, penalty = penaltyb1,
        min = wb1, bqicu_min = bqicu1[wb1], vars = best_sub[[wb1]]
    )
    res_cs_bqicu[[j]] <- list(
        bqicu = bqicu2, gof = gof2, penalty = penaltyb2,
        min = wb2, bqicu_min = bqicu2[wb2], vars = best_sub[[wb2]]
    )
    res_ar1_bqicu[[j]] <- list(
        bqicu = bqicu3, gof = gof3, penalty = penaltyb3,
        min = wb3, bqicu_min = bqicu3[wb3], vars = best_sub[[wb3]]
    )
    setTxtProgressBar(pb, j)
    if (j == nsims) {
        close(pb)
    }
}

# Save simulations ####
save(
    res_indp_qic, res_cs_qic, res_ar1_qic,
    res_indp_qicu, res_cs_qicu, res_ar1_qicu,
    res_indp_bqicu, res_cs_bqicu, res_ar1_bqicu,
    file = "./simulations/outputs/norm_glsx_bqicu_sim/norm_glsx_bqicu_sim.RData"
)
