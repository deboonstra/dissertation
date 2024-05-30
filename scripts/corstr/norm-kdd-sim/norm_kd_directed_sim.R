# The function of this script file is to run a new simulation investigating
# the selection properties and distribution of Kullback-Leibler directed
# divergence. This simulation is an extension of the simulations ran in
# norm_modified_kdd_sim, as we are using the model-based covariance matrix
# based on the working correlation structure.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Setup functions
invert <- if ("MASS" %in% loadedNamespaces()) {
  MASS::ginv
} else {
  solve
}

# Creating output sub-directory ####
sub_dir <- "./outputs/corstr/norm-kdd-sim/norm-kd-directed-sim/"
if (!dir.exists(sub_dir)) {
  dir.create(sub_dir)
}

# Defining global data simulation settings ####
tol <- .Machine$double.eps
nsims <- 1000L
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
form <- stats::as.formula(
  paste0("y~", paste0("X", which(beta != 0), collapse = "+"))
)
l <- sum(beta != 0) + 1 # kb limit
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

# Simulation for kb values ####
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
      # Simulating data ####
      dat <- sim_data(
        N = res_basis$N[ell], n = res_basis$n[ell], beta = beta, rho = rho,
        corstr = res_basis$corstr[ell]
      )

      ## Converting to data.frame ####
      dat <- data.frame(y = dat$y, dat$X, id = dat$id)

      # Fitting models and getting kb values ####
      kb <- rep(NA, times = length(work_corstr)) ## statistics
      kb0 <- NA ## oracle
      names(kb) <- names_work_corstr
      for (k in seq_along(work_corstr)) {
        ## Fitting model ####
        fit <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[k],
          Mv = mv[k]
        )

        ### Getting initial kb values ####
        kb[k] <- kd(object = fit)

        ##### Obtaining new kb values if kb < 0 and kb > l * 10
        iter <- 0
        while ((is.na(kb[k]) || kb[k] < 0 || kb[k] > l * 10) && iter <= 100) {
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

          # Getting kb ####
          kb[k] <- kd(object = fit)

          # Updating iteration count ####
          iter <- iter + 1
        }

        ## Getting final kb values ####
        ## A kb value will be missing if there are still negative variances
        ## or too large variances
        ### Model-based orcale values ####
        if ((is.na(kb[k]) || kb[k] < 0 || kb[k] > l * 10) && iter == 100) {
          kb[k] <- NA
        } else {
          kb[k] <- kb[k]
        }
      }

      ## Obtaining the oracle kb ####
      ## based on the minimum kb statistic
      if (all(!is.na(kb))) {
        fit0 <- gee::gee(
          formula = form, id = id,
          data = dat,
          corstr = work_corstr[which.min(kb)],
          Mv = mv[which.min(kb)],
          maxiter = 100
        )

        # Calculating components of divergence ####

        ##  Number mean parameter estimates ####
        k <- length(fit0$coefficients)

        ## Covariance matrices ####
        mb <- fit0$naive.variance

        ## Calculating trace terms ####
        tr <- sum(diag(crossprod(x = invert(mb, tol = tol), y = sigma0[[ell]])))

        ## log-determinant terms ####
        ld <- log(det(mb) / det(sigma0[[ell]]))

        # Calculating KL divergence
        kb0 <- 0.5 * (ld + tr - k)
      } else {
        kb0 <- NA
      }

      # Comparison of KL divergence values ####
      # There will be four comparisons, which are outlined below.
      # 1. Exchangeable and AR(1)
      # 2. Independence, exchangeable, and AR(1)
      # 3. Independence, exchangeable, AR(1), and AR(3)
      # 4. Independence, exchangeable, AR(1), AR(3), and unstructured
      # If a kb value is missing for a comparison, then NO correlation
      # structure will be selected. A missing value will be reported as a fair
      # comparison can not be made.

      ## 1. Exchangeable and AR(1) only ####
      include_corstr <- c("exchangeable", "AR-M")
      kb_sel1 <- kb[which(work_corstr %in% include_corstr & mv == 1)]
      sel1 <- ifelse(
        test = sum(is.na(kb_sel1)) == 0,
        yes = names(kb_sel1)[which.min(kb_sel1)],
        no = NA
      )

      ## 2. Independence, exchangeable, and AR(1) ####
      include_corstr <- c("independence", "exchangeable", "AR-M")
      kb_sel2 <- kb[which(work_corstr %in% include_corstr & mv == 1)]
      sel2 <- ifelse(
        test = sum(is.na(kb_sel2)) == 0,
        yes = names(kb_sel2)[which.min(kb_sel2)],
        no = NA
      )

      ## 3. Independence, exchangeable, AR(1), and AR(3) ####
      kb_sel3 <- kb[which(work_corstr != "unstructured")]
      sel3 <- ifelse(
        test = sum(is.na(kb_sel3)) == 0,
        yes = names(kb_sel3)[which.min(kb_sel3)],
        no = NA
      )

      ## 4. Independence, exchangeable, AR(1), AR(3), unstructured ####
      sel4 <- ifelse(
        test = sum(is.na(kb)) == 0,
        yes = names(kb)[which.min(kb)],
        no = NA
      )

      # Creating output ####
      out <- data.frame(
        sims = j,
        N = res_basis$N[ell],
        n = res_basis$n[ell],
        corstr = res_basis$corstr[ell],
        matrix(
          data = kb,
          nrow = 1, ncol = length(work_corstr),
          dimnames = list(NULL, names_work_corstr), byrow = TRUE
        ),
        kb0 = kb0,
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
  file = paste0(sub_dir, "norm_kd_directed_sim.rds")
)