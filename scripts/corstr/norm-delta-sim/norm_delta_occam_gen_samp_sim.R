# The function of this script file is to determine the selection properties
# of a new correlation structure selection measure called "delta", where
# Occam's window is implemented. The window is determined by randomly
# generating data and obtaining the information criterion sampling distribution
# for the largest correlation structure of interest.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
sub_dir <- "./outputs/corstr/norm-delta-sim/norm-delta-occam-gen-samp-sim/"
if (!dir.exists(sub_dir)) {
  dir.create(sub_dir)
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
form <- stats::as.formula(
  paste0("y~", paste0("X", which(beta != 0), collapse = "+"))
)
l <- sum(beta != 0) + 1 # delta limit
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
strict_upper <- FALSE

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
cat("Obtaining the true covariance matrix\n")
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

# Obtaining the sampling distributions
cat("Obtaining the sampling distributions\n")
set.seed(1997)
dc <- parallel::detectCores()
cl <- parallel::makeCluster(spec = (dc - 1))
un_delta <- vector(mode = "list", length = nrow(res_basis))
names(un_delta) <- paste0(
  "N", res_basis$N, "n", res_basis$n, "corstr_", res_basis$corstr
)
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
  un_delta[[ell]] <- parallel::parSapply(
    cl = cl,
    X = seq_len(nsims),
    FUN = function(j) {
      # Simulating data ####
      dat <- occam_data(
        N = res_basis$N[ell], n = res_basis$n[ell], p = l
      )
      ## Converting to data.frame ####
      dat <- data.frame(y = dat$y, dat$X, id = dat$id)
      # Fitting models and getting delta values ####

      ## Fitting model ####
      fit <- gee::gee(
        formula = form, id = id,
        data = dat,
        corstr = "unstructured",
        Mv = 1
      )
      ### Pulling covariance matrices ####
      mb <- fit$naive.variance
      vr <- fit$robust.variance

      ### Getting initial delta values ####
      dd <- delta(a = mb, b = vr)

      ##### Obtaining new delta values if delta < l and delta > l * 10
      iter <- 0
      while ((is.na(dd) || dd < l || dd > l * 10) && iter <= 100) {
        # Simulating data ####
        dat <- occam_data(
          N = res_basis$N[ell], n = res_basis$n[ell], p = l
        )
        ## Converting to data.frame ####
        dat <- data.frame(y = dat$y, dat$X, id = dat$id)

        # Fitting model ####
        fit <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = "unstructured",
          Mv = 1
        )

        ## Pulling covariance matrices ####
        mb <- fit$naive.variance
        vr <- fit$robust.variance

        # Getting delta ####
        dd <- delta(a = mb, b = vr)

        # Updating iteration count ####
        iter <- iter + 1
      }
      if ((is.na(dd) || dd < l || dd > l * 10) && iter == 100) {
        dd <- NA
      } else {
        dd <- dd
      }
      # Returning output ####
      return(dd)
    }
  )
  # Updating progress ####
  cat("\n")
}

# Simulation for delta values ####
cat("Simulation for delta values\n")
set.seed(1997)
res <- vector(mode = "list", length = nrow(res_basis))
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
      # Simulating data ####
      dat <- sim_data(
        N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
        corstr = res_basis$corstr[ell]
      )

      ## Converting to data.frame ####
      dat <- data.frame(y = dat$y, dat$X, id = dat$id)

      # Fitting models and getting delta values ####
      dd <- rep(NA, times = length(work_corstr)) ## statistics
      dd0 <- NA ## oracle
      names(dd) <- names_work_corstr
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

        ### Getting initial delta values ####
        dd[k] <- delta(a = mb, b = vr)

        ##### Obtaining new delta values if delta < l and delta > l * 10
        iter <- 0
        while ((is.na(dd[k]) || dd[k] < l || dd[k] > l * 10) && iter <= 100) {
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

          ## Pulling covariance matrices ####
          mb <- fit$naive.variance
          vr <- fit$robust.variance

          # Getting delta ####
          dd[k] <- delta(a = mb, b = vr)

          # Updating iteration count ####
          iter <- iter + 1
        }

        ## Getting final delta values ####
        ## A delta value will be missing if there are still negative variances
        ## or too large variances
        ### Model-based orcale values ####
        if ((is.na(dd[k]) || dd[k] < l || dd[k] > l * 10) && iter == 100) {
          dd[k] <- NA
        } else {
          dd[k] <- dd[k]
        }
      }

      ## Obtaining the oracle delta ####
      ## based on the minimum delta statistic
      if (all(!is.na(dd))) {
        fit0 <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[which.min(dd)],
          Mv = mv[which.min(dd)],
          maxiter = 100
        )
        dd0 <- delta(a = fit0$naive.variance, b = sigma0[[ell]])
      } else {
        dd0 <- NA
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
      ## Because exchangeable and AR(1) have the same number of parameters, we
      ## don't have to implement Occam's window.
      include_corstr <- c("exchangeable", "AR-M")
      delta_sel1 <- dd[which(work_corstr %in% include_corstr & mv == 1)]
      mv_sel1 <- mv[which(work_corstr %in% include_corstr & mv == 1)]
      sel1 <- ifelse(
        test = sum(is.na(delta_sel1)) == 0,
        yes = names(delta_sel1)[which.min(delta_sel1)],
        no = NA
      )

      ## 2. Independence, exchangeable, and AR(1) ####
      ## The window here is determined by the generating correlation structure,
      ## as in practice you would be able to see if the correlation among
      ## observations taper off (i.e., AR(1)) or constant (i.e., exchangeable)
      include_corstr <- c("independence", "exchangeable", "AR-M")
      delta_sel2 <- dd[which(work_corstr %in% include_corstr & mv == 1)]
      names(delta_sel2) <- work_corstr[which(work_corstr %in% include_corstr & mv == 1)] #nolint
      mv_sel2 <- mv[which((work_corstr %in% include_corstr) & (mv == 1))]
      if (sum(is.na(delta_sel2)) == 0) {
        occam_sel2 <- occam(
          x = delta_sel2,
          ideal = l,
          window = unname(quantile(x = un_delta[[ell]], probs = 0.9)) - l,
          mv = mv_sel2,
          strict_upper = strict_upper
        )
        sel2 <- occam_sel2$name
        sel2 <- ifelse(
          test = sel2 == "AR-M",
          yes = paste0("ar", occam_sel2$nparms),
          no = sel2
        )
      } else {
        sel2 <- NA
      }

      ## 3. Independence, exchangeable, AR(1), and AR(3) ####
      delta_sel3 <- dd[which(work_corstr != "unstructured")]
      names(delta_sel3) <- work_corstr[which(work_corstr != "unstructured")] #nolint
      mv_sel3 <- mv[which(work_corstr != "unstructured")]
      if (sum(is.na(delta_sel3)) == 0) {
        occam_sel3 <- occam(
          x = delta_sel3,
          ideal = l,
          window = unname(quantile(x = un_delta[[ell]], probs = 0.9)) - l,
          mv = mv_sel3,
          strict_upper = strict_upper
        )
        sel3 <- occam_sel3$name
        sel3 <- ifelse(
          test = sel3 == "AR-M",
          yes = paste0("ar", occam_sel3$nparms),
          no = sel3
        )
      } else {
        sel3 <- NA
      }

      ## 4. Independence, exchangeable, AR(1), AR(3), unstructured ####
      if (sum(is.na(dd)) == 0) {
        names(dd) <- work_corstr
        occam_sel4 <- occam(
          x = dd,
          ideal = l,
          window = unname(quantile(x = un_delta[[ell]], probs = 0.9)) - l,
          mv = mv,
          strict_upper = strict_upper
        )
        sel4 <- occam_sel4$name
        sel4 <- ifelse(
          test = sel4 == "AR-M",
          yes = paste0("ar", occam_sel3$nparms),
          no = sel4
        )
      } else {
        sel4 <- NA
      }

      # Creating output ####
      out <- data.frame(
        sims = j,
        N = res_basis$N[ell],
        n = res_basis$n[ell],
        corstr = res_basis$corstr[ell],
        matrix(
          data = dd,
          nrow = 1, ncol = length(work_corstr),
          dimnames = list(NULL, names_work_corstr), byrow = TRUE
        ),
        delta0 = dd0,
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
  file = paste0(sub_dir, "norm_delta_occam_gen_samp_sim.rds")
)