# The function of this script is to run a simulation investigating the selection
# properties of KDD penalized by expected optimism (eo). Additionally, we will
# be determining how to properly specify eo in these simulations. The mean
# difference component will not be calculated in these simulations as the
# limiting distributions of all the models will have the same mean vector.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output directory if one does not exists
if (!dir.exists("./outputs/corstr/norm-kdd-sim/norm-kdd-eo-sim/")) {
  dir.create("./outputs/corstr/norm-kdd-sim/norm-kdd-eo-sim/")
}

# Defining global data simulation settings ####
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
nsims <- 1000
wcorstr <- c("exchangeable", "ar1", "independence", "unstructured")

# Simulation for original KDD -- KDD(I) ####
# What this means is that the model-based covariance matrix used in KDD is based
# on the independence working correlation structure similar to QIC and CIC.
set.seed(1997)
res <- data.frame(
  wcorstr = rep(x = wcorstr, each = nsims),
  kdd_unpenalized = NA,
  kdd_y = NA,
  kdd1 = NA,
  eo1 = NA,
  kdd_final1 = NA,
  kdd2 = NA,
  eo2 = NA,
  kdd_final2 = NA
)
dc <- parallel::detectCores()
cl <- parallel::makeCluster(dc - 1)
pb <- utils::txtProgressBar(min = 0, max = nrow(res), style = 3)
for (j in seq_len(nrow(res))) {
  # Data
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = "exchangeable")
  y <- dat$y
  X <- dat$X
  id <- dat$id
  dat <- data.frame(y = y, X, id = id)

  # Fit model
  fit <- geepack::geeglm(
    y ~ X1 + X2 + X3,
    id = id, data = dat, corstr = res$wcorstr[j]
  )

  # Fit smallest model
  mf <- geepack::geeglm(
    y ~ 1,
    id = id, data = dat, corstr = "independence"
  )

  # eo ####
  r <- 100 # replicates
  dat_ystar <- fit$data
  dat_zstar <- fit$data
  parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
  ## Calculating replicate expected optimism values ####
  eo <- parallel::parLapply(
    cl = cl,
    X = seq_len(r),
    fun = function(i) {
      dat_ystar$y <- rnorm(n = length(mf$y), mean = coef(mf))
      dat_zstar$y <- rnorm(n = length(mf$y), mean = coef(mf))
      mod_y <- geepack::geeglm(
        formula = fit$formula, data = dat_ystar, id = id, corstr = fit$corstr
      )
      mod_z <- geepack::geeglm(
        formula = fit$formula, data = dat_zstar, id = id, corstr = fit$corstr
      )
      kdd_y <- kdd(object1 = mod_y, object2 = mod_y)
      kdd1 <- kdd(object1 = mod_y, object2 = mod_z)
      kdd2 <- kdd(object1 = mod_z, object2 = mod_y)
      eo1 <- kdd1 - kdd_y
      eo2 <- kdd2 - kdd_y
      out <- data.frame(
        kdd_y = kdd_y,
        kdd1 = kdd1,
        eo1 = eo1,
        kdd2 = kdd2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final KDD values with expected optimism penalty
  res$kdd_unpenalized[j] <- kdd(object1 = fit, object2 = fit)
  res$kdd_final1[j] <- res$kdd_unpenalized[j] + res$eo1[j]
  res$kdd_final2[j] <- res$kdd_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_kdd <- res

# Simulation for modified KDD -- KDD(R) ####
# What this means is that the model-based covariance matrix used in KDD is based
# on the working correlation structure similar to modified CIC.
set.seed(1997)
res <- data.frame(
  wcorstr = rep(x = wcorstr, each = nsims),
  kdd_unpenalized = NA,
  kdd_y = NA,
  kdd1 = NA,
  eo1 = NA,
  kdd_final1 = NA,
  kdd2 = NA,
  eo2 = NA,
  kdd_final2 = NA
)
dc <- parallel::detectCores()
cl <- parallel::makeCluster(dc - 1)
pb <- utils::txtProgressBar(min = 0, max = nrow(res), style = 3)
for (j in seq_len(nrow(res))) {
  # Data
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = "exchangeable")
  y <- dat$y
  X <- dat$X
  id <- dat$id
  dat <- data.frame(y = y, X, id = id)

  # Fit model
  fit <- geepack::geeglm(
    y ~ X1 + X2 + X3,
    id = id, data = dat, corstr = res$wcorstr[j]
  )

  # Fit smallest model
  mf <- geepack::geeglm(
    y ~ 1,
    id = id, data = dat, corstr = "independence"
  )

  # eo ####
  r <- 100 # replicates
  dat_ystar <- fit$data
  dat_zstar <- fit$data
  parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
  ## Calculating replicate expected optimism values ####
  eo <- parallel::parLapply(
    cl = cl,
    X = seq_len(r),
    fun = function(i) {
      dat_ystar$y <- rnorm(n = length(mf$y), mean = coef(mf))
      dat_zstar$y <- rnorm(n = length(mf$y), mean = coef(mf))
      mod_y <- geepack::geeglm(
        formula = fit$formula, data = dat_ystar, id = id, corstr = fit$corstr
      )
      mod_z <- geepack::geeglm(
        formula = fit$formula, data = dat_zstar, id = id, corstr = fit$corstr
      )
      kdd_y <- kdd(object1 = mod_y, object2 = mod_y, modified = TRUE)
      kdd1 <- kdd(object1 = mod_y, object2 = mod_z, modified = TRUE)
      kdd2 <- kdd(object1 = mod_z, object2 = mod_y, modified = TRUE)
      eo1 <- kdd1 - kdd_y
      eo2 <- kdd2 - kdd_y
      out <- data.frame(
        kdd_y = kdd_y,
        kdd1 = kdd1,
        eo1 = eo1,
        kdd2 = kdd2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final KDD values with expected optimism penalty
  res$kdd_unpenalized[j] <- kdd(object1 = fit, object2 = fit, modified = TRUE)
  res$kdd_final1[j] <- res$kdd_unpenalized[j] + res$eo1[j]
  res$kdd_final2[j] <- res$kdd_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_modified_kdd <- res

# Exporting ####
save(
  res_kdd, res_modified_kdd,
  file = "./outputs/corstr/norm-kdd-sim/norm-kdd-eo-sim/sim_res.RData"
)