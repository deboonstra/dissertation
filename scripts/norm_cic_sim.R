# The function of script file is run a simulation to investigate the selection
# properties of CIC. For information about this simulation, see the meeting
# notes for 2023-09-07 and the document explaining this simulation:
# `outputs/norm_cic_sim/norm_cic_sim.R`.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 1000L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- c("exchangeable", "ar1")
work_corstr <- c("independence", "exchangeable", "ar1", "unstructured")

# Simulation ####
set.seed(1997)
res_cic1_cs <- vector(mode = "list", length = nsims)
res_cic1_ar1 <- vector(mode = "list", length = nsims)
res_cic2_cs <- vector(mode = "list", length = nsims)
res_cic2_ar1 <- vector(mode = "list", length = nsims)
res_cic3_cs <- vector(mode = "list", length = nsims)
res_cic3_ar1 <- vector(mode = "list", length = nsims)
pb <- utils::txtProgressBar(min = 0, max = nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  ## Two data sets are being simulated for different correlation structures
  ## to generate the random component
  ### Exchangeable
  dat_cs <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr[1])
  ### AR(1)
  dat_ar1 <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr[2])
  ### Converting to data frames
  dat_cs <- data.frame(y = dat_cs$y, dat_cs$X, id = dat_cs$id)
  dat_ar1 <- data.frame(y = dat_ar1$y, dat_ar1$X, id = dat_ar1$id)

  ## Fitting model and getting CIC values ####
  cic_cs <- rep(NA, length(work_corstr)) # for exchangeable data
  cic_ar1 <- rep(NA, length(work_corstr)) # for AR(1) data
  names(cic_cs) <- names(cic_ar1) <- work_corstr
  for (k in seq_along(work_corstr)) {
    ### Fitting models
    fit_cs <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_cs, corstr = work_corstr[k]
    )
    fit_ar1 <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_ar1, corstr = work_corstr[k]
    )
    ### Getting CIC values
    cic_cs[k] <- cic(fit_cs)
    cic_ar1[k] <- cic(fit_ar1)
  }

  ## Comparison  of CIC values ####
  ## There will be three comparisons, which are outlined below.
  ## 1. Exchangeable and AR(1)
  ## 2. Independence, exchangeable, and AR(1)
  ## 3. Independence, exchangeable, AR(1), and unstructured

  ### 1. Exchangeable and AR(1) only ####
  cic1_cs <- cic_cs[which(work_corstr %in% c("exchangeable", "ar1"))]
  cic1_ar1 <- cic_ar1[which(work_corstr %in% c("exchangeable", "ar1"))]
  w1_cs <- which.min(cic1_cs)
  w1_ar1 <- which.min(cic1_ar1)

  ### 2. Independence, exchangeable, and AR(1)
  cic2_cs <- cic_cs[which(work_corstr != "unstructured")]
  cic2_ar1 <- cic_ar1[which(work_corstr != "unstructured")]
  w2_cs <- which.min(cic2_cs)
  w2_ar1 <- which.min(cic2_ar1)

  ### 3. Independence, exchangeable, AR(1), and unstructured
  cic3_cs <- cic_cs
  cic3_ar1 <- cic_ar1
  w3_cs <- which.min(cic3_cs)
  w3_ar1 <- which.min(cic3_ar1)

  ## Updating result ####
  res_cic1_cs[[j]] <- list(
    ic = cic1_cs, corstr = names(cic1_cs),
    ic_min = cic1_cs[w1_cs], corstr_min = names(cic1_cs)[w1_cs]
  )
  res_cic1_ar1[[j]] <- list(
    ic = cic1_ar1, corstr = names(cic1_ar1),
    ic_min = cic1_ar1[w1_ar1], corstr_min = names(cic1_ar1)[w1_ar1]
  )
  res_cic2_cs[[j]] <- list(
    ic = cic2_cs, corstr = names(cic2_cs),
    ic_min = cic2_cs[w2_cs], corstr_min = names(cic2_cs)[w2_cs]
  )
  res_cic2_ar1[[j]] <- list(
    ic = cic2_ar1, corstr = names(cic2_ar1),
    ic_min = cic2_ar1[w2_ar1], corstr_min = names(cic2_ar1)[w2_ar1]
  )
  res_cic3_cs[[j]] <- list(
    ic = cic3_cs, corstr = names(cic3_cs),
    ic_min = cic3_cs[w1_cs], corstr_min = names(cic3_cs)[w3_cs]
  )
  res_cic3_ar1[[j]] <- list(
    ic = cic3_ar1, corstr = names(cic3_ar1),
    ic_min = cic3_ar1[w3_ar1], corstr_min = names(cic3_ar1)[w3_ar1]
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation results ####
if (!dir.exists("./outputs/norm-cic-sim")) {
  dir.create("./outputs/norm-cic-sim")
}

save(
  res_cic1_cs, res_cic1_ar1,
  res_cic2_cs, res_cic2_ar1,
  res_cic3_cs, res_cic3_ar1,
  file = "./outputs/norm-cic-sim/norm_cic_sim.RData"
)