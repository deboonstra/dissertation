# The function of script file is run a simulation to investigate the selection
# properties of the modified inverted logarithm of the determinants. For
# information about this simulation, see the meeting notes for 2023-10-05 and
# the document explaining this simulation:
# `outputs/norm-modified-inv-log-det-sim/norm_modified_inv_log_det_sim.Rmd`.

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
res_ld1_cs <- vector(mode = "list", length = nsims)
res_ld1_ar1 <- vector(mode = "list", length = nsims)
res_ld2_cs <- vector(mode = "list", length = nsims)
res_ld2_ar1 <- vector(mode = "list", length = nsims)
res_ld3_cs <- vector(mode = "list", length = nsims)
res_ld3_ar1 <- vector(mode = "list", length = nsims)
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

  ## Fitting model and getting ld values ####
  ld_cs <- rep(NA, length(work_corstr)) # for exchangeable data
  ld_ar1 <- rep(NA, length(work_corstr)) # for AR(1) data
  names(ld_cs) <- names(ld_ar1) <- work_corstr
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
    ### Getting ld values
    ld_cs[k] <- -1 * log_det(fit_cs, modified = TRUE)
    ld_ar1[k] <- -1 * log_det(fit_ar1, modified = TRUE)
  }

  ## Comparison  of ld values ####
  ## There will be three comparisons, which are outlined below.
  ## 1. Exchangeable and AR(1)
  ## 2. Independence, exchangeable, and AR(1)
  ## 3. Independence, exchangeable, AR(1), and unstructured

  ### 1. Exchangeable and AR(1) only ####
  ld1_cs <- ld_cs[which(work_corstr %in% c("exchangeable", "ar1"))]
  ld1_ar1 <- ld_ar1[which(work_corstr %in% c("exchangeable", "ar1"))]
  w1_cs <- which.min(ld1_cs)
  w1_ar1 <- which.min(ld1_ar1)

  ### 2. Independence, exchangeable, and AR(1)
  ld2_cs <- ld_cs[which(work_corstr != "unstructured")]
  ld2_ar1 <- ld_ar1[which(work_corstr != "unstructured")]
  w2_cs <- which.min(ld2_cs)
  w2_ar1 <- which.min(ld2_ar1)

  ### 3. Independence, exchangeable, AR(1), and unstructured
  ld3_cs <- ld_cs
  ld3_ar1 <- ld_ar1
  w3_cs <- which.min(ld3_cs)
  w3_ar1 <- which.min(ld3_ar1)

  ## Updating result ####
  res_ld1_cs[[j]] <- list(
    ic = ld1_cs, corstr = names(ld1_cs),
    ic_min = ld1_cs[w1_cs], corstr_min = names(ld1_cs)[w1_cs]
  )
  res_ld1_ar1[[j]] <- list(
    ic = ld1_ar1, corstr = names(ld1_ar1),
    ic_min = ld1_ar1[w1_ar1], corstr_min = names(ld1_ar1)[w1_ar1]
  )
  res_ld2_cs[[j]] <- list(
    ic = ld2_cs, corstr = names(ld2_cs),
    ic_min = ld2_cs[w2_cs], corstr_min = names(ld2_cs)[w2_cs]
  )
  res_ld2_ar1[[j]] <- list(
    ic = ld2_ar1, corstr = names(ld2_ar1),
    ic_min = ld2_ar1[w2_ar1], corstr_min = names(ld2_ar1)[w2_ar1]
  )
  res_ld3_cs[[j]] <- list(
    ic = ld3_cs, corstr = names(ld3_cs),
    ic_min = ld3_cs[w1_cs], corstr_min = names(ld3_cs)[w3_cs]
  )
  res_ld3_ar1[[j]] <- list(
    ic = ld3_ar1, corstr = names(ld3_ar1),
    ic_min = ld3_ar1[w3_ar1], corstr_min = names(ld3_ar1)[w3_ar1]
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation results ####
if (!dir.exists("./outputs/corstr/norm-ksd-sim/norm-modified-inv-log-det-sim")) { #nolint
  dir.create("./outputs/corstr/norm-ksd-sim/norm-modified-inv-log-det-sim")
}

save(
  res_ld1_cs, res_ld1_ar1,
  res_ld2_cs, res_ld2_ar1,
  res_ld3_cs, res_ld3_ar1,
  file = paste0(
    "./outputs/corstr/norm-ksd-sim/norm-modified-inv-log-det-sim/",
    "norm_modified_inv_log_det_sim.RData"
  )
)
