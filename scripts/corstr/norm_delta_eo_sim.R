# The function of this script file is to determine the selection properties
# of a new correlation structure selection measure called "delta" that includes
# a penalty based on expected optimism.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta-eo-sim/")) {
  dir.create("./outputs/corstr/norm-delta-eo-sim/")
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- 200
n <- c(3, 4, 5, 6, 7, 8, 9, 10)
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
  cat(paste0("Simulations for n = ", n[ell], "\n"))
  pb <- utils::txtProgressBar(
    min = 0, max = nsims_cov, style = 3, initial = "\n"
  )
  for (i in seq_len(nsims_cov)) {
    ## Simulating data ####
    ## Two data sets are being simulated for different correlation structures
    ## to generate the random component
    ### Exchangeable
    dat_cs <- sim_data(
      N = N, n = n[ell], beta = beta, rho = rho, corstr = corstr[1]
    )
    ### AR(1)
    dat_ar1 <- sim_data(
      N = N, n = n[ell], beta = beta, rho = rho, corstr = corstr[2]
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
names(sigma0_cs) <- names(sigma0_ar1) <- paste0("n", n)

# Simulation ####
set.seed(1997)
res_delta1_cs <- vector(mode = "list", length = length(n) * nsims)
res_delta1_ar1 <- vector(mode = "list", length = length(n) * nsims)
res_delta2_cs <- vector(mode = "list", length = length(n) * nsims)
res_delta2_ar1 <- vector(mode = "list", length = length(n) * nsims)
res_delta3_cs <- vector(mode = "list", length = length(n) * nsims)
res_delta3_ar1 <- vector(mode = "list", length = length(n) * nsims)
res_delta0_cs <- vector(mode = "list", length = length(n))
res_delta0_ar1 <- vector(mode = "list", length = length(n))
res_eo_cs <- vector(mode = "list", length = length(n))
res_eo_ar1 <- vector(mode = "list", length = length(n))
dc <- parallel::detectCores()
cl <- parallel::makeCluster(dc - 1)
for (ell in seq_along(n)) {
  res_delta0_cs[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  res_delta0_ar1[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  res_eo_cs[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  res_eo_ar1[[ell]] <- matrix(
    data = NA, nrow = nsims, ncol = length(work_corstr),
    dimnames = list(paste0("Simulation: ", seq_len(nsims)), work_corstr)
  )
  cat(paste0("Simulations for n = ", n[ell], "\n"))
  pb <- utils::txtProgressBar(min = 0, max = nsims, style = 3)
  for (j in seq_len(nsims)) {
    ## Simulating data ####
    ## Two data sets are being simulated for different correlation structures
    ## to generate the random component
    ### Exchangeable
    dat_cs <- sim_data(
      N = N, n = n[ell], beta = beta, rho = rho, corstr = corstr[1]
    )
    ### AR(1)
    dat_ar1 <- sim_data(
      N = N, n = n[ell], beta = beta, rho = rho, corstr = corstr[2]
    )
    ### Converting to data frames
    dat_cs <- data.frame(y = dat_cs$y, dat_cs$X, id = dat_cs$id)
    dat_ar1 <- data.frame(y = dat_ar1$y, dat_ar1$X, id = dat_ar1$id)

    ## Fitting model and getting delta values ####
    delta_cs <- rep(NA, length(work_corstr)) # for exchangeable data
    delta_ar1 <- rep(NA, length(work_corstr)) # for AR(1) data
    names(delta_cs) <- names(delta_ar1) <- work_corstr
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

      ### Expected optimism calculations ####

      #### Number of replications ####
      r <- 100

      #### Generating replications based on the smallest model ####
      yy_cs <- matrix(
        data = rnorm(
          n = r * length(fit_cs$y),
          mean = 0,
          sd = 1
        ),
        nrow = r, ncol = length(fit_cs$y)
      )
      zz_cs <- matrix(
        data = rnorm(
          n = r * length(fit_cs$y),
          mean = 0,
          sd = 1
        ),
        nrow = r, ncol = length(fit_cs$y)
      )

      yy_ar1 <- matrix(
        data = rnorm(
          n = r * length(fit_ar1$y),
          mean = 0,
          sd = 1
        ),
        nrow = r, ncol = length(fit_ar1$y)
      )
      zz_ar1 <- matrix(
        data = rnorm(
          n = r * length(fit_ar1$y),
          mean = 0,
          sd = 1
        ),
        nrow = r, ncol = length(fit_ar1$y)
      )

      #### Cluster export ####
      parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))

      #### eo calculations ####
      eo <- parallel::parLapply(
        cl = cl,
        X = seq_len(r),
        fun = function(x) {
          # Replication model fits ####
          mod_y_cs <- geepack::geeglm(
            formula = yy_cs[x, ] ~ 1, id = fit_cs$id, corstr = fit_cs$corstr
          )
          mod_z_cs <- geepack::geeglm(
            formula = zz_cs[x, ] ~ 1, id = fit_cs$id, corstr = fit_cs$corstr
          )
          mod_y_ar1 <- geepack::geeglm(
            formula = yy_ar1[x, ] ~ 1, id = fit_ar1$id, corstr = fit_ar1$corstr
          )
          mod_z_ar1 <- geepack::geeglm(
            formula = zz_ar1[x, ] ~ 1, id = fit_ar1$id, corstr = fit_ar1$corstr
          )
          # Calculation of delta values ####
          delta_yy_cs <- delta(
            a = mod_y_cs$geese$vbeta.naiv,
            b = mod_y_cs$geese$vbeta
          )
          delta_yz_cs <- delta(
            a = mod_y_cs$geese$vbeta.naiv,
            b = mod_z_cs$geese$vbeta
          )
          delta_yy_ar1 <- delta(
            a = mod_y_ar1$geese$vbeta.naiv,
            b = mod_y_ar1$geese$vbeta
          )
          delta_yz_ar1 <- delta(
            a = mod_y_ar1$geese$vbeta.naiv,
            b = mod_z_ar1$geese$vbeta
          )
          # Calculation of eo ####
          eo_cs <- delta_yz_cs - delta_yy_cs
          eo_ar1 <- delta_yz_ar1 - delta_yy_ar1

          # Return ####
          return(data.frame(eo_cs = eo_cs, eo_ar1 = eo_ar1))
        }
      )
      eo <- colMeans(dplyr::bind_rows(eo))

      ### Getting delta values
      delta_cs[k] <- delta(
        a = fit_cs$geese$vbeta.naiv,
        b = fit_cs$geese$vbeta,
        penalty = eo[1]
      )
      delta_ar1[k] <- delta(
        a = fit_ar1$geese$vbeta.naiv,
        b = fit_ar1$geese$vbeta,
        penalty = eo[2]
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

      #### Getting eo values ####

      ##### Exchangeable
      res_eo_cs[[ell]][j, k] <- eo[1]

      ##### AR(1)
      res_delta0_ar1[[ell]][j, k] <- eo[2]
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
      n = n[ell], ic = delta1_cs, corstr = names(delta1_cs),
      ic_min = delta1_cs[w1_cs], corstr_min = names(delta1_cs)[w1_cs]
    )
    res_delta1_ar1[[j]] <- list(
      n = n[ell], ic = delta1_ar1, corstr = names(delta1_ar1),
      ic_min = delta1_ar1[w1_ar1], corstr_min = names(delta1_ar1)[w1_ar1]
    )
    res_delta2_cs[[j]] <- list(
      n = n[ell], ic = delta2_cs, corstr = names(delta2_cs),
      ic_min = delta2_cs[w2_cs], corstr_min = names(delta2_cs)[w2_cs]
    )
    res_delta2_ar1[[j]] <- list(
      n = n[ell], ic = delta2_ar1, corstr = names(delta2_ar1),
      ic_min = delta2_ar1[w2_ar1], corstr_min = names(delta2_ar1)[w2_ar1]
    )
    res_delta3_cs[[j]] <- list(
      n = n[ell], ic = delta3_cs, corstr = names(delta3_cs),
      ic_min = delta3_cs[w1_cs], corstr_min = names(delta3_cs)[w3_cs]
    )
    res_delta3_ar1[[j]] <- list(
      n = n[ell], ic = delta3_ar1, corstr = names(delta3_ar1),
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
parallel::stopCluster(cl)
names(res_delta0_cs) <- names(res_delta0_ar1) <- paste0("n", n)
names(res_eo_cs) <- names(res_eo_ar1) <- paste0("n", n)
res_names <- paste0(
  "n",
  rep(n, each = nsims),
  "_",
  rep(seq_len(nsims), times = length(n))
)
names(res_delta1_cs) <- names(res_delta2_ar1) <- res_names
names(res_delta2_cs) <- names(res_delta2_ar1) <- res_names
names(res_delta3_cs) <- names(res_delta3_ar1) <- res_names

# Save simulation results ####
save(
  res_delta1_cs, res_delta1_ar1,
  res_delta2_cs, res_delta2_ar1,
  res_delta3_cs, res_delta3_ar1,
  res_delta0_cs, res_delta0_ar1,
  res_eo_cs, res_eo_ar1,
  file = "./outputs/corstr/norm-delta-eo-sim/norm_delta_eo_sim.RData"
)