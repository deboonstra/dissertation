# The function of script file is run a simulation to investigate the selection
# properties of KDD. For information about this simulation, see the meeting
# notes for 2023-09-21 and the document explaining this simulation:
# `outputs/norm_kdd_sim/norm_kdd_sim.Rmd`.

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
res_kl1_cs <- vector(mode = "list", length = nsims)
res_kl1_ar1 <- vector(mode = "list", length = nsims)
res_kl2_cs <- vector(mode = "list", length = nsims)
res_kl2_ar1 <- vector(mode = "list", length = nsims)
res_kl3_cs <- vector(mode = "list", length = nsims)
res_kl3_ar1 <- vector(mode = "list", length = nsims)
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

  ## Fitting model and getting kl values ####
  kl_cs <- rep(NA, length(work_corstr)) # for exchangeable data
  kl_ar1 <- rep(NA, length(work_corstr)) # for AR(1) data
  names(kl_cs) <- names(kl_ar1) <- work_corstr
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
    ### Getting kl values
    kl_cs[k] <- kl(fit_cs)
    kl_ar1[k] <- kl(fit_ar1)
  }

  ## Comparison  of kl values ####
  ## There will be three comparisons, which are outlined below.
  ## 1. Exchangeable and AR(1)
  ## 2. Independence, exchangeable, and AR(1)
  ## 3. Independence, exchangeable, AR(1), and unstructured

  ### 1. Exchangeable and AR(1) only ####
  kl1_cs <- kl_cs[which(work_corstr %in% c("exchangeable", "ar1"))]
  kl1_ar1 <- kl_ar1[which(work_corstr %in% c("exchangeable", "ar1"))]
  w1_cs <- which.min(kl1_cs)
  w1_ar1 <- which.min(kl1_ar1)

  ### 2. Independence, exchangeable, and AR(1)
  kl2_cs <- kl_cs[which(work_corstr != "unstructured")]
  kl2_ar1 <- kl_ar1[which(work_corstr != "unstructured")]
  w2_cs <- which.min(kl2_cs)
  w2_ar1 <- which.min(kl2_ar1)

  ### 3. Independence, exchangeable, AR(1), and unstructured
  kl3_cs <- kl_cs
  kl3_ar1 <- kl_ar1
  w3_cs <- which.min(kl3_cs)
  w3_ar1 <- which.min(kl3_ar1)

  ## Updating result ####
  res_kl1_cs[[j]] <- list(
    ic = kl1_cs, corstr = names(kl1_cs),
    ic_min = kl1_cs[w1_cs], corstr_min = names(kl1_cs)[w1_cs]
  )
  res_kl1_ar1[[j]] <- list(
    ic = kl1_ar1, corstr = names(kl1_ar1),
    ic_min = kl1_ar1[w1_ar1], corstr_min = names(kl1_ar1)[w1_ar1]
  )
  res_kl2_cs[[j]] <- list(
    ic = kl2_cs, corstr = names(kl2_cs),
    ic_min = kl2_cs[w2_cs], corstr_min = names(kl2_cs)[w2_cs]
  )
  res_kl2_ar1[[j]] <- list(
    ic = kl2_ar1, corstr = names(kl2_ar1),
    ic_min = kl2_ar1[w2_ar1], corstr_min = names(kl2_ar1)[w2_ar1]
  )
  res_kl3_cs[[j]] <- list(
    ic = kl3_cs, corstr = names(kl3_cs),
    ic_min = kl3_cs[w1_cs], corstr_min = names(kl3_cs)[w3_cs]
  )
  res_kl3_ar1[[j]] <- list(
    ic = kl3_ar1, corstr = names(kl3_ar1),
    ic_min = kl3_ar1[w3_ar1], corstr_min = names(kl3_ar1)[w3_ar1]
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation results ####
if (!dir.exists("./outputs/norm-kdd-sim")) {
  dir.create("./outputs/norm-kdd-sim")
}

save(
  res_kl1_cs, res_kl1_ar1,
  res_kl2_cs, res_kl2_ar1,
  res_kl3_cs, res_kl3_ar1,
  file = "./outputs/norm-kdd-sim/norm_kdd_sim.RData"
)
