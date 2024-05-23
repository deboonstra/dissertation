# The function of script file is run a simulation to investigate the selection
# properties of CIC,
# where the beta coefficients are large.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
sub_dir <- paste0(
  "./outputs/corstr/norm-cic-sim/",
  "norm-cic-large-beta-sim/"
)
if (!dir.exists(sub_dir)) {
  dir.create(sub_dir)
}
# Defining global data simulation settings ####
nsims <- 1000L
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
beta <- c(2000, 3000, 500, 0, 0, 0)
form <- stats::as.formula(
  paste0("y~", paste0("X", which(beta != 0), collapse = "+"))
)
l <- sum(beta != 0) + 1 # cic limit
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

# Simulation for cic values ####
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

      # Fitting models and getting cic values ####
      cc <- rep(NA, times = length(work_corstr))
      names(cc) <- names_work_corstr
      for (k in seq_along(work_corstr)) {
        ## Fitting model ####
        fit <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[k],
          Mv = mv[k]
        )

        ### Getting initial cic values ####
        cc[k] <- gee::cic(object = fit)

        ##### Obtaining new cic values if CIC < 0 and CIC > l * 10
        iter <- 0
        while ((is.na(cc[k]) || cc[k] < 0 || cc[k] > l * 10) && iter <= 100) {
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

          # Getting cic ####
          cc[k] <- gee::cic(object = fit)

          # Updating iteration count ####
          iter <- iter + 1
        }

        ## Getting final cic values ####
        ## A cic value will be missing if there are still negative variances
        ## or too large variances
        ### Model-based orcale values ####
        if ((is.na(cc[k]) || cc[k] < 0 || cc[k] > l * 10) && iter == 100) {
          cc[k] <- NA
        } else {
          cc[k] <- cc[k]
        }
      }

      ## Obtaining the oracle CIC ####
      ## based on the minimum cic statistic
      if (all(!is.na(cc))) {
        # Re-fitting best model ####
        fit0 <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[which.min(cc)],
          Mv = mv[which.min(cc)],
          maxiter = 100
        )

        # Obtaining covariance matrices ####

        ### Re-fitting independence model
        model_indep <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = "independence",
          Mv = 1,
          maxiter = 100,
          b = fit0$coefficients,
          scale.fix = TRUE,
          scale.value = fit0$scale
        )

        ## Getting covariance matrices ####

        ### Inverse of model-based covariance under independence ####
        omega <- MASS::ginv(
          X = model_indep$naive.variance, tol = .Machine$double.eps
        )

        ### Robust covariance given working correlation structure ####
        vr <- fit0$robust.variance

        # Calculating CIC ####
        cc0 <- sum(diag(x = crossprod(x = omega, y = vr)))
      } else {
        cc0 <- NA
      }

      # Comparison of cic values ####
      # There will be four comparisons, which are outlined below.
      # 1. Exchangeable and AR(1)
      # 2. Independence, exchangeable, and AR(1)
      # 3. Independence, exchangeable, AR(1), and AR(3)
      # 4. Independence, exchangeable, AR(1), AR(3), and unstructured
      # If a cic value is missing for a comparison, then NO correlation
      # structure will be selected. A missing value will be reported as a fair
      # comparison can not be made.

      ## 1. Exchangeable and AR(1) only ####
      include_corstr <- c("exchangeable", "AR-M")
      cic_sel1 <- cc[which(work_corstr %in% include_corstr & mv == 1)]
      sel1 <- ifelse(
        test = sum(is.na(cic_sel1)) == 0,
        yes = names(cic_sel1)[which.min(cic_sel1)],
        no = NA
      )

      ## 2. Independence, exchangeable, and AR(1) ####
      include_corstr <- c("independence", "exchangeable", "AR-M")
      cic_sel2 <- cc[which(work_corstr %in% include_corstr & mv == 1)]
      sel2 <- ifelse(
        test = sum(is.na(cic_sel2)) == 0,
        yes = names(cic_sel2)[which.min(cic_sel2)],
        no = NA
      )

      ## 3. Independence, exchangeable, AR(1), and AR(3) ####
      cic_sel3 <- cc[which(work_corstr != "unstructured")]
      sel3 <- ifelse(
        test = sum(is.na(cic_sel3)) == 0,
        yes = names(cic_sel3)[which.min(cic_sel3)],
        no = NA
      )

      ## 4. Independence, exchangeable, AR(1), AR(3), unstructured ####
      sel4 <- ifelse(
        test = sum(is.na(cc)) == 0,
        yes = names(cc)[which.min(cc)],
        no = NA
      )

      # Creating output ####
      out <- data.frame(
        sims = j,
        N = res_basis$N[ell],
        n = res_basis$n[ell],
        corstr = res_basis$corstr[ell],
        matrix(
          data = cc,
          nrow = 1, ncol = length(work_corstr),
          dimnames = list(NULL, names_work_corstr), byrow = TRUE
        ),
        cic0 = cc0,
        sel1 = sel1,
        sel2 = sel2,
        sel3 = sel3,
        sel4 = sel4
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
  file = paste0(sub_dir, "norm_cic_large_beta_sim.rds")
)