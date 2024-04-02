# The function of this script file is to run a simulation investigating the
# distribution and selection properties of the oracle forms of "delta". This
# is a deeper dive than what was done in the `norm_delta_sim` simulations.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta0-sim/")) {
  dir.create("./outputs/corstr/norm-delta0-sim/")
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
form <- stats::as.formula(
  paste0("y~", paste0("X", which(beta != 0), collapse = "+"))
)
rho <- 0.5
corstr <- c("exchangeable", "ar1")
work_corstr <- c(
  "independence",
  "exchangeable",
  "AR-M",
  "AR-M",
  "unstructured"
)
names_work_corstr <- c(
  "independence",
  "exchangeable",
  "ar1",
  "ar3",
  "unstructured"
)
mv <- c(1, 1, 1, 3, 1)

# Creating a basis data.frame to store the simulations results ####
res_basis <- data.frame(
  sims = rep(x = nsims, each = length(corstr)),
  N = rep(x = N, each = length(corstr)),
  n = rep(x = n, each = length(corstr)),
  corstr = rep(x = corstr, times = length(N)),
  matrix(
    data = NA,
    nrow = length(corstr) * length(N), ncol = length(work_corstr),
    dimnames = list(NULL, names_work_corstr)
  )
)

# Obtaining the true covariance matrix ####
set.seed(1997)
nsims_cov <- nsims * 10L
D <- 0 # nolint: Reducing nsims_cov if diag(mb) < 0
sigma0 <- vector(mode = "list", length = nrow(res_basis))
names(sigma0) <- paste0(
  "N", res_basis$N, "n", res_basis$n, "corstr_", res_basis$corstr
)
for (ell in seq_along(sigma0)) {
  sigma0[[ell]] <- matrix(
    data = 0,
    nrow = length(which(beta != 0)) + 1,
    ncol = length(which(beta != 0)) + 1
  )
  cat(
    paste0(
      "Simulations for N = ",
      res_basis$N[ell],
      " and n = ",
      res_basis$n[ell],
      " given corstr = ",
      res_basis$corstr[ell],
      "\n"
    )
  )
  pb <- utils::txtProgressBar(
    min = 0, max = nsims_cov, style = 3, initial = "\n"
  )
  for (i in seq_len(nsims_cov)) {
    # Simulating data ####
    dat <- sim_data(
      N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
      corstr = res_basis$corstr[ell]
    )
    ## Converting to data.frame ####
    dat <- data.frame(y = dat$y, dat$X, id = dat$id)

    # Fitting model ####
    fit <- gee::gee(
      formula = form, id = id,
      data = dat,
      corstr = ifelse(
        test = res_basis$corstr[ell] == "ar1",
        yes = "AR-M",
        no = res_basis$corstr[ell]
      )
    )

    ## Pulling model-based covariance matrix ####
    mb <- fit$naive.variance

    ### Obtaining new covariance matrice if diag(mb) < 0 ####
    iter <- 0
    while (sum(diag(mb) < 0) >= 1 && iter <= 100) {
      # Simulating data ####
      dat <- sim_data(
        N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
        corstr = res_basis$corstr[ell]
      )
      ## Converting to data.frame ####
      dat <- data.frame(y = dat$y, dat$X, id = dat$id)

      # Fitting model ####
      fit <- gee::gee(
        formula = form, id = id,
        data = dat,
        corstr = ifelse(
          test = res_basis$corstr[ell] == "ar1",
          yes = "AR-M",
          no = res_basis$corstr[ell]
        )
      )

      ## Pulling model-based covariance matrix ####
      mb <- fit$naive.variance

      # Updating iteration count ####
      iter <- iter + 1
    }

    # Updating reduction variable ####
    if (iter == 100 && sum(diag(mb) < 0) >= 1) {
      diag(mb) <- matrix(data = 0, nrow = dim(mb)[1], nrol = dim(mb)[2])
      D <- D + 1 # nolint
    }

    # Summing model-based covariance matrices ####
    sigma0[[ell]] <- sigma0[[ell]] + mb

    # Updating progress bar ####
    utils::setTxtProgressBar(pb, i)
    if (i == nsims_cov) {
      close(pb)
      cat("\n")
    }
  }
  # Averaging model-based covariance matrices ####
  sigma0[[ell]] <- (1 / (nsims_cov - D)) * sigma0[[ell]]
}

# Simulation for orcale delta values ####
set.seed(1997)
res <- vector(mode = "list", length = nrow(res_basis))
dc <- parallel::detectCores()
cl <- parallel::makeCluster(spec = (dc - 1))
for (ell in seq_len(nrow(res_basis))) {
  cat(
    paste0(
      "Simulations for N = ",
      res_basis$N[ell],
      " and n = ",
      res_basis$n[ell],
      " given corstr = ",
      res_basis$corstr[ell],
      "\n"
    )
  )
  parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
  res[[ell]] <- parallel::parLapply(
    cl = cl,
    X = seq_len(res_basis$sims[ell]),
    fun = function(j) {
      # Simulation data ####
      dat <- sim_data(
        N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
        corstr = res_basis$corstr[ell]
      )

      ## Converting to data.frame ####
      dat <- data.frame(y = dat$y, dat$X, id = dat$id)

      # Fitting models and getting delta values ####
      delta0_mb <- delta0_vr <- rep(NA, times = length(work_corstr))
      names(delta0_mb) <- names(delta0_vr) <- names_work_corstr
      for (k in seq_along(work_corstr)) {
        ## Fitting model ####
        fit <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[k],
          Mv = mv[k]
        )

        ### Pulling covariance matrices ####
        mb <- fit$naive.variance
        vr <- fit$robust.variance

        #### Obtaining new covariance matrice if diag(mb) < 0 ####
        iter_mb <- 0
        while (sum(diag(mb) < 0) >= 1 && iter_mb <= 100) {
          # Simulating data ####
          dat <- sim_data(
            N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
            corstr = res_basis$corstr[ell]
          )
          ## Converting to data.frame ####
          dat <- data.frame(y = dat$y, dat$X, id = dat$id)

          # Fitting model ####
          fit <- gee::gee(
            formula = form, id = id,
            data = dat,
            corstr = work_corstr[k],
            Mv = mv[k]
          )

          ## Pulling model-based covariance matrix ####
          mb <- fit$naive.variance

          # Updating iteration count ####
          iter_mb <- iter_mb + 1
        }

        #### Obtaining new covariance matrice if diag(vr) < 0 ####
        iter_vr <- 0
        while (sum(diag(vr) < 0) >= 1 && iter_vr <= 100) {
          # Simulating data ####
          dat <- sim_data(
            N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
            corstr = res_basis$corstr[ell]
          )
          ## Converting to data.frame ####
          dat <- data.frame(y = dat$y, dat$X, id = dat$id)

          # Fitting model ####
          fit <- gee::gee(
            formula = form, id = id,
            data = dat,
            corstr = work_corstr[k],
            Mv = mv[k]
          )

          ## Pulling model-based covariance matrix ####
          vr <- fit$robust.variance

          # Updating iteration count ####
          iter_vr <- iter_vr + 1
        }

        ## Getting delta values ####
        ## A delta value will be missing if there are still negative variances
        ### Model-based orcale values ####
        if (iter_mb == 100 && sum(diag(mb) < 0) >= 1) {
          delta0_mb[k] <- NA
        } else {
          delta0_mb[k] <- delta(
            a = mb,
            b = sigma0[[ell]]
          )
        }

        ### Robust orcale values ####
        if (iter_vr == 100 && sum(diag(vr) < 0) >= 1) {
          delta0_vr[k] <- NA
        } else {
          delta0_vr[k] <- delta(
            a = vr,
            b = sigma0[[ell]]
          )
        }
      }

      # Comparison of delta values ####
      # There will be four comparisons, which are outlined below.
      # 1. Exchangeable and AR(1)
      # 2. Independence, exchangeable, and AR(1)
      # 3. Independence, exchangeable, AR(1), and AR(3)
      # 4. Independence, exchangeable, AR(1), AR(3), and unstructured
      # If a delta value is missing for a comparison, then NO correlation
      # structure will be selected. A missing value will be reported as a fair
      # comparison can not be made.

      ## 1. Exchangeable and AR(1) only ####
      include_corstr <- c("exchangeable", "AR-M")
      delta0_mb1 <- delta0_mb[which(work_corstr %in% include_corstr & mv == 1)]
      sel_mb1 <- ifelse(
        test = sum(is.na(delta0_mb1)) == 0,
        yes = names(delta0_mb1)[which.min(delta0_mb1)],
        no = NA
      )
      delta0_vr1 <- delta0_vr[which(work_corstr %in% include_corstr & mv == 1)]
      sel_vr1 <- ifelse(
        test = sum(is.na(delta0_vr1)) == 0,
        yes = names(delta0_vr1)[which.min(delta0_vr1)],
        no = NA
      )

      ## 2. Independence, exchangeable, and AR(1) ####
      include_corstr <- c("independence", "exchangeable", "AR-M")
      delta0_mb2 <- delta0_mb[which(work_corstr %in% include_corstr & mv == 1)]
      sel_mb2 <- ifelse(
        test = sum(is.na(delta0_mb2)) == 0,
        yes = names(delta0_mb2)[which.min(delta0_mb2)],
        no = NA
      )
      delta0_vr2 <- delta0_vr[which(work_corstr %in% include_corstr & mv == 1)]
      sel_vr2 <- ifelse(
        test = sum(is.na(delta0_vr2)) == 0,
        yes = names(delta0_vr2)[which.min(delta0_vr2)],
        no = NA
      )

      ## 3. Independence, exchangeable, AR(1), and AR(3) ####
      delta0_mb3 <- delta0_mb[which(work_corstr != "unstructured")]
      sel_mb3 <- ifelse(
        test = sum(is.na(delta0_mb3)) == 0,
        yes = names(delta0_mb3)[which.min(delta0_mb3)],
        no = NA
      )
      delta0_vr3 <- delta0_vr[which(work_corstr != "unstructured")]
      sel_vr3 <- ifelse(
        test = sum(is.na(delta0_vr3)) == 0,
        yes = names(delta0_vr3)[which.min(delta0_vr3)],
        no = NA
      )

      ## 4. Independence, exchangeable, AR(1), AR(3), unstructured ####
      sel_mb4 <- ifelse(
        test = sum(is.na(delta0_mb)) == 0,
        yes = names(delta0_mb)[which.min(delta0_mb)],
        no = NA
      )
      sel_vr4 <- ifelse(
        test = sum(is.na(delta0_vr)) == 0,
        yes = names(delta0_vr)[which.min(delta0_vr)],
        no = NA
      )

      # Creating output ####
      out <- data.frame(
        sims = rep(x = j, times = 2),
        N = rep(x = res_basis$N[ell], times = 2),
        n = rep(x = res_basis$n[ell], times = 2),
        corstr = rep(x = res_basis$corstr[ell], times = 2),
        type = c("mb", "vr"),
        matrix(
          data = c(delta0_mb, delta0_vr),
          nrow = 2, ncol = length(work_corstr),
          dimnames = list(NULL, names_work_corstr), byrow = TRUE
        ),
        sel1 = c(sel_mb1, sel_vr1),
        sel2 = c(sel_mb2, sel_vr2),
        sel3 = c(sel_mb3, sel_vr3),
        sel4 = c(sel_mb4, sel_vr4)
      )

      # Returning output ####
      return(out)
    }
  )
  res[[ell]] <- dplyr::bind_rows(res[[ell]])
  # Updating progress ####
  cat("\n")
}

# Shutting down clusters ####
parallel::stopCluster(cl = cl)

# Converting list of simulation results to a data.frame of results ####
res <- dplyr::bind_rows(res)

# Exporting simulation results ####
saveRDS(
  object = res,
  file = "./outputs/corstr/norm-delta0-sim/norm_delta0_sim.rds"
)