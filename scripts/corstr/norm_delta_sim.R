# The function of this script file is to determine the selection properties
# of a new correlation structure selection measure called "delta".

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta-sim/")) {
  dir.create("./outputs/corstr/norm-delta-sim/")
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- c(rep(200, 8), 160, 100, 80, 50, 40, 32)
n <- c(3, 4, 5, 6, 7, 8, 9, 10, 5, 8, 10, 16, 20, 25)
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- c("exchangeable", "ar1")
work_corstr <- c("independence", "exchangeable", "ar1", "unstructured")

# Obtaining true covariance matrix ####
set.seed(1997)
nsims_cov <- nsims * 10L
sigma0_cs <- vector(mode = "list", length = length(n))
sigma0_ar1 <- vector(mode = "list", length = length(n))
for (ell in seq_along(n)) {
  sigma0_cs[[ell]] <- matrix(
    data = 0,
    nrow = length(which(beta != 0)) + 1,
    ncol = length(which(beta != 0)) + 1
  )
  sigma0_ar1[[ell]] <- matrix(
    data = 0,
    nrow = length(which(beta != 0)) + 1,
    ncol = length(which(beta != 0)) + 1
  )
  cat(paste0("Simulations for N = ", N[ell], " and n = ", n[ell], "\n"))
  pb <- utils::txtProgressBar(
    min = 0, max = nsims_cov, style = 3, initial = "\n"
  )
  for (i in seq_len(nsims_cov)) {
    ## Simulating data ####
    ## Two data sets are being simulated for different correlation structures
    ## to generate the random component
    ### Exchangeable
    dat_cs <- sim_data(
      N = N[ell], n = n[ell], beta = beta, rho = rho, corstr = corstr[1]
    )
    ### AR(1)
    dat_ar1 <- sim_data(
      N = N[ell], n = n[ell], beta = beta, rho = rho, corstr = corstr[2]
    )
    ### Converting to data frames
    dat_cs <- data.frame(y = dat_cs$y, dat_cs$X, id = dat_cs$id)
    dat_ar1 <- data.frame(y = dat_ar1$y, dat_ar1$X, id = dat_ar1$id)

    ## Fitting models ####
    ### Exchangeable
    fit_cs <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_cs, corstr = corstr[1]
    )
    mb_cs <- fit_cs$geese$vbeta.naiv
    ### AR(1)
    fit_ar1 <- geepack::geeglm(
      formula = y ~ X1 + X2 + X3, id = id,
      data = dat_ar1, corstr = corstr[2]
    )
    mb_ar1 <- fit_ar1$geese$vbeta.naiv

    ## Pulling model-based covariance matrices ####
    sigma0_cs[[ell]] <- sigma0_cs[[ell]] + ((1 / nsims_cov) * mb_cs)
    sigma0_ar1[[ell]] <- sigma0_ar1[[ell]] + ((1 / nsims_cov) * mb_ar1)
    ## Updating progress bar ####
    utils::setTxtProgressBar(pb, i)
    if (i == nsims_cov) {
      close(pb)
      cat("\n")
    }
  }
}
names(sigma0_cs) <- names(sigma0_ar1) <- paste0("N", N, "n", n)


# Simulation ####
set.seed(1997)
res_delta1_cs <- vector(mode = "list", length = length(n))
res_delta1_ar1 <- vector(mode = "list", length = length(n))
res_delta2_cs <- vector(mode = "list", length = length(n))
res_delta2_ar1 <- vector(mode = "list", length = length(n))
res_delta3_cs <- vector(mode = "list", length = length(n))
res_delta3_ar1 <- vector(mode = "list", length = length(n))
res_delta0_cs <- vector(mode = "list", length = length(n))
res_delta0_ar1 <- vector(mode = "list", length = length(n))
for (ell in seq_along(n)) {
  res_delta0_cs[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  res_delta0_ar1[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  res_delta1_cs[[ell]] <- vector(mode = "list", length = nsims)
  res_delta1_ar1[[ell]] <- vector(mode = "list", length = nsims)
  res_delta2_cs[[ell]] <- vector(mode = "list", length = nsims)
  res_delta2_ar1[[ell]] <- vector(mode = "list", length = nsims)
  res_delta3_cs[[ell]] <- vector(mode = "list", length = nsims)
  res_delta3_ar1[[ell]] <- vector(mode = "list", length = nsims)
  cat(paste0("Simulations for N = ", N[ell], " and n = ", n[ell], "\n"))
  pb <- utils::txtProgressBar(
    min = 0, max = nsims, style = 3, initial = "\n"
  )
  for (j in seq_len(nsims)) {
    ## Simulating data ####
    ## Two data sets are being simulated for different correlation structures
    ## to generate the random component
    ### Exchangeable
    dat_cs <- sim_data(
      N = N[ell], n = n[ell], beta = beta, rho = rho, corstr = corstr[1]
    )
    ### AR(1)
    dat_ar1 <- sim_data(
      N = N[ell], n = n[ell], beta = beta, rho = rho, corstr = corstr[2]
    )
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
      fit_ar1 <- geepack::geeglm(
        formula = y ~ X1 + X2 + X3, id = id,
        data = dat_ar1, corstr = work_corstr[k]
      )
      ### Getting delta values ####
      delta_cs[k] <- delta(
        a = fit_cs$geese$vbeta.naiv,
        b = fit_cs$geese$vbeta
      )
      delta_ar1[k] <- delta(
        a = fit_ar1$geese$vbeta.naiv,
        b = fit_ar1$geese$vbeta
      )

      #### Getting true delta values ####

      ##### Exchangeable
      res_delta0_cs[[ell]][j, k] <- delta(
        a = fit_cs$geese$vbeta.naiv,
        b = sigma0_cs[[ell]]
      )

      ##### AR(1)
      res_delta0_ar1[[ell]][j, k] <- delta(
        a = fit_ar1$geese$vbeta.naiv,
        b = sigma0_ar1[[ell]]
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
    res_delta1_cs[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta1_cs, corstr = names(delta1_cs),
      ic_min = delta1_cs[w1_cs], corstr_min = names(delta1_cs)[w1_cs]
    )
    res_delta1_ar1[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta1_ar1, corstr = names(delta1_ar1),
      ic_min = delta1_ar1[w1_ar1], corstr_min = names(delta1_ar1)[w1_ar1]
    )
    res_delta2_cs[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta2_cs, corstr = names(delta2_cs),
      ic_min = delta2_cs[w2_cs], corstr_min = names(delta2_cs)[w2_cs]
    )
    res_delta2_ar1[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta2_ar1, corstr = names(delta2_ar1),
      ic_min = delta2_ar1[w2_ar1], corstr_min = names(delta2_ar1)[w2_ar1]
    )
    res_delta3_cs[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta3_cs, corstr = names(delta3_cs),
      ic_min = delta3_cs[w3_cs], corstr_min = names(delta3_cs)[w3_cs]
    )
    res_delta3_ar1[[ell]][[j]] <- list(
      N = N[ell], n = n[ell], ic = delta3_ar1, corstr = names(delta3_ar1),
      ic_min = delta3_ar1[w3_ar1], corstr_min = names(delta3_ar1)[w3_ar1]
    )
    ## Updating progress bar ####
    utils::setTxtProgressBar(pb, j)
    if (j == nsims) {
      close(pb)
      cat("\n")
    }
  }
}
names(res_delta0_cs) <- names(res_delta0_ar1) <- paste0("N", N, "n", n)
names(res_delta1_cs) <- names(res_delta2_ar1) <- paste0("N", N, "n", n)
names(res_delta2_cs) <- names(res_delta2_ar1) <- paste0("N", N, "n", n)
names(res_delta3_cs) <- names(res_delta3_ar1) <- paste0("N", N, "n", n)

# Save simulation results ####
save(
  res_delta1_cs, res_delta1_ar1,
  res_delta2_cs, res_delta2_ar1,
  res_delta3_cs, res_delta3_ar1,
  res_delta0_cs, res_delta0_ar1,
  file = "./outputs/corstr/norm-delta-sim/norm_delta_sim.RData"
)