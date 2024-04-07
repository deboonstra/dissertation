# The function of this script file is to determine the selection properties
# of a new correlation structure selection measure called "delta" that includes
# a penalty for having multiple correlation parameters to estimate.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta-sim/norm-delta-penalty-sim/")) {
  dir.create("./outputs/corstr/norm-delta-sim/norm-delta-penalty-sim/")
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- c("exchangeable", "ar1")
work_corstr <- c("independence", "exchangeable", "ar1", "unstructured")

# Obtaining true covariance matrix ####
set.seed(1997)
nsims_cov <- nsims * 10L
sigma0_cs <- matrix(
  data = 0,
  nrow = length(which(beta != 0)) + 1,
  ncol = length(which(beta != 0)) + 1
)
sigma0_ar1 <- matrix(
  data = 0,
  nrow = length(which(beta != 0)) + 1,
  ncol = length(which(beta != 0)) + 1
)
pb <- utils::txtProgressBar(min = 0, max = nsims_cov, style = 3)
for (i in seq_len(nsims_cov)) {
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

  ## Fitting models ####
  ### Exchangeable
  fit_cs <- geepack::geeglm(
    formula = y ~ X1 + X2 + X3, id = id,
    data = dat_cs, corstr = corstr[1]
  )
  ### AR(1)
  fit_ar1 <- geepack::geeglm(
    formula = y ~ X1 + X2 + X3, id = id,
    data = dat_ar1, corstr = corstr[2]
  )

  ## Pulling model-based covariance matrices ####
  sigma0_cs <- sigma0_cs + ((1 / nsims_cov) * fit_cs$geese$vbeta.naiv)
  sigma0_ar1 <- sigma0_ar1 + ((1 / nsims_cov) * fit_ar1$geese$vbeta.naiv)
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, i)
  if (i == nsims_cov) {
    close(pb)
    cat("\n")
  }
}

# Simulation ####
set.seed(1997)
res_delta1_cs <- vector(mode = "list", length = nsims)
res_delta1_ar1 <- vector(mode = "list", length = nsims)
res_delta2_cs <- vector(mode = "list", length = nsims)
res_delta2_ar1 <- vector(mode = "list", length = nsims)
res_delta3_cs <- vector(mode = "list", length = nsims)
res_delta3_ar1 <- vector(mode = "list", length = nsims)
res_delta0_cs <- matrix(
  data = NA, nrow = nsims, ncol = length(work_corstr),
  dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
)
res_delta0_ar1 <- matrix(
  data = NA, nrow = nsims, ncol = length(work_corstr),
  dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
)
pb <- utils::txtProgressBar(min = 0, max = nsims, initial = "\n", style = 3)
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

  ## Fitting model and getting delta values ####
  delta_cs <- rep(NA, length(work_corstr)) # for exchangeable data
  delta_ar1 <- rep(NA, length(work_corstr)) # for AR(1) data
  names(delta_cs) <- names(delta_ar1) <- work_corstr
  for (k in seq_along(work_corstr)) {
    ### Fitting models ####
    fit_cs <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_cs, corstr = work_corstr[k]
    )
    penalty_cs <- ifelse(
      test = length(fit_cs$geese$alpha) == 0,
      yes = 0,
      no = -1 / length(fit_cs$geese$alpha)
    )
    fit_ar1 <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_ar1, corstr = work_corstr[k]
    )
    penalty_ar1 <- ifelse(
      test = length(fit_ar1$geese$alpha) == 0,
      yes = 0,
      no = -1 / length(fit_ar1$geese$alpha)
    )
    ### Getting delta values ####
    delta_cs[k] <- delta(
      a = fit_cs$geese$vbeta.naiv,
      b = fit_cs$geese$vbeta,
      penalty = penalty_cs
    )
    delta_ar1[k] <- delta(
      a = fit_ar1$geese$vbeta.naiv,
      b = fit_ar1$geese$vbeta,
      penalty = penalty_ar1
    )

    #### Getting true delta values ####

    ##### Exchangeable
    res_delta0_cs[j, k] <- delta(
      a = fit_cs$geese$vbeta.naiv,
      b = sigma0_cs,
      penalty = penalty_cs
    )

    ##### AR(1)
    res_delta0_ar1[j, k] <- delta(
      a = fit_ar1$geese$vbeta.naiv,
      b = sigma0_ar1,
      penalty = penalty_ar1
    )
  }

  ## Comparison  of delta values ####
  ## There will be three comparisons, which are outlined below.
  ## 1. Exchangeable and AR(1)
  ## 2. Independence, exchangeable, and AR(1)
  ## 3. Independence, exchangeable, AR(1), and unstructured

  ### 1. Exchangeable and AR(1) only ####
  delta1_cs <- delta_cs[which(work_corstr %in% c("exchangeable", "ar1"))]
  delta1_ar1 <- delta_ar1[which(work_corstr %in% c("exchangeable", "ar1"))]
  w1_cs <- which.min(delta1_cs)
  w1_ar1 <- which.min(delta1_ar1)

  ### 2. Independence, exchangeable, and AR(1)
  delta2_cs <- delta_cs[which(work_corstr != "unstructured")]
  delta2_ar1 <- delta_ar1[which(work_corstr != "unstructured")]
  w2_cs <- which.min(delta2_cs)
  w2_ar1 <- which.min(delta2_ar1)

  ### 3. Independence, exchangeable, AR(1), and unstructured
  delta3_cs <- delta_cs
  delta3_ar1 <- delta_ar1
  w3_cs <- which.min(delta3_cs)
  w3_ar1 <- which.min(delta3_ar1)

  ## Updating result ####
  res_delta1_cs[[j]] <- list(
    ic = delta1_cs, corstr = names(delta1_cs),
    ic_min = delta1_cs[w1_cs], corstr_min = names(delta1_cs)[w1_cs]
  )
  res_delta1_ar1[[j]] <- list(
    ic = delta1_ar1, corstr = names(delta1_ar1),
    ic_min = delta1_ar1[w1_ar1], corstr_min = names(delta1_ar1)[w1_ar1]
  )
  res_delta2_cs[[j]] <- list(
    ic = delta2_cs, corstr = names(delta2_cs),
    ic_min = delta2_cs[w2_cs], corstr_min = names(delta2_cs)[w2_cs]
  )
  res_delta2_ar1[[j]] <- list(
    ic = delta2_ar1, corstr = names(delta2_ar1),
    ic_min = delta2_ar1[w2_ar1], corstr_min = names(delta2_ar1)[w2_ar1]
  )
  res_delta3_cs[[j]] <- list(
    ic = delta3_cs, corstr = names(delta3_cs),
    ic_min = delta3_cs[w1_cs], corstr_min = names(delta3_cs)[w3_cs]
  )
  res_delta3_ar1[[j]] <- list(
    ic = delta3_ar1, corstr = names(delta3_ar1),
    ic_min = delta3_ar1[w3_ar1], corstr_min = names(delta3_ar1)[w3_ar1]
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
    cat("\n")
  }
}

# Save simulation results ####
save(
  res_delta1_cs, res_delta1_ar1,
  res_delta2_cs, res_delta2_ar1,
  res_delta3_cs, res_delta3_ar1,
  res_delta0_cs, res_delta0_ar1,
  file = paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-penalty-sim/norm_delta_penalty_sim.RData"
  )
)