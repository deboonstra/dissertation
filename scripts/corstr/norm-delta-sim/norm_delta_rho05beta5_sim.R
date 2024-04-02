# The function of this script file is to determine the selection properties
# of a new correlation structure selection measure called "delta".
# However, there will be five non-zero mean parameters in these simulations and
# the correlation coefficient will be 0.05.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta-sim/norm-delta-rho05beta5-sim/")) {
  dir.create("./outputs/corstr/norm-delta-sim/norm-delta-rho05beta5-sim/")
}

# Defining global data simulation settings ####
nsims <- 1000L
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
beta <- c(2.0, 3.0, 0.5, 1.0, -1.0, 0)
form <- stats::as.formula(
  paste0("y~", paste0("X", which(beta != 0), collapse = "+"))
)
rho <- 0.05
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
      dd <- rep(NA, times = length(work_corstr))
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
        cond_mb <- (iter_mb == 100 && sum(diag(mb) < 0) >= 1)
        cond_vr <- (iter_vr == 100 && sum(diag(vr) < 0) >= 1)
        if (cond_mb || cond_vr) {
          dd[k] <- NA
        } else {
          dd[k] <- delta(a = mb, b = vr)
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
      delta_sel1 <- dd[which(work_corstr %in% include_corstr & mv == 1)]
      sel1 <- ifelse(
        test = sum(is.na(delta_sel1)) == 0,
        yes = names(delta_sel1)[which.min(delta_sel1)],
        no = NA
      )

      ## 2. Independence, exchangeable, and AR(1) ####
      include_corstr <- c("independence", "exchangeable", "AR-M")
      delta_sel2 <- dd[which(work_corstr %in% include_corstr & mv == 1)]
      sel2 <- ifelse(
        test = sum(is.na(delta_sel2)) == 0,
        yes = names(delta_sel2)[which.min(delta_sel2)],
        no = NA
      )

      ## 3. Independence, exchangeable, AR(1), and AR(3) ####
      delta_sel3 <- dd[which(work_corstr != "unstructured")]
      sel3 <- ifelse(
        test = sum(is.na(delta_sel3)) == 0,
        yes = names(delta_sel3)[which.min(delta_sel3)],
        no = NA
      )

      ## 4. Independence, exchangeable, AR(1), AR(3), unstructured ####
      sel4 <- ifelse(
        test = sum(is.na(dd)) == 0,
        yes = names(dd)[which.min(dd)],
        no = NA
      )

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
  file = paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-rho05beta5-sim/norm_delta_rho05beta5_sim.rds"
  )
)