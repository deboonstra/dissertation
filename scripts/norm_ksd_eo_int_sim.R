# The function of this script is to run a simulation investigating the selection
# properties of KSD penalized by expected optimism (eo). Additionally, we will
# be determining how to properly specify eo in these simulations. The mean
# difference component will not be calculated in these simulations as the
# limiting distributions of all the models will have the same mean vector.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output directory if one does not exists
if (!dir.exists("./outputs/norm-ksd-eo-int-sim/")) {
  dir.create("./outputs/norm-ksd-eo-int-sim/")
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
  ksd_unpenalized = NA,
  ksd_y = NA,
  ksd1 = NA,
  eo1 = NA,
  ksd_final1 = NA,
  ksd2 = NA,
  eo2 = NA,
  ksd_final2 = NA
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

  # eo ####
  r <- 100 # replicates

  ## Generating replicates based on the smallest model ####
  ## Smallest model is an intercept only model with an independence working
  ## correlation structure. Thus, we generate any iid N(mu, sigma), where
  ## mu and sigma are constant and fixed.
  yy <- matrix(
    data = rnorm(
      n = r * length(fit$y),
      mean = 0,
      sd = 1
    ),
    nrow = r, ncol = length(fit$y)
  )
  zz <- matrix(
    data = rnorm(
      n = r * length(fit$y),
      mean = 0,
      sd = 1
    ),
    nrow = r, ncol = length(fit$y)
  )
  parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
  ## Calculating replicate expected optimism values ####
  eo <- parallel::parLapply(
    cl = cl,
    X = seq_len(r),
    fun = function(x) {
      mod_y <- geepack::geeglm(
        formula = yy[x, ] ~ 1, id = fit$id, corstr = fit$corstr
      )
      mod_z <- geepack::geeglm(
        formula = zz[x, ] ~ 1, id = fit$id, corstr = fit$corstr
      )
      ksd_y <- kdd(
        object1 = mod_y, object2 = mod_y,
        modified = FALSE, symmetric = TRUE
      )
      ksd1 <- kdd(
        object1 = mod_y, object2 = mod_z,
        modified = FALSE, symmetric = TRUE
      )
      ksd2 <- kdd(
        object1 = mod_z, object2 = mod_y,
        modified = FALSE, symmetric = TRUE
      )
      eo1 <- ksd1 - ksd_y
      eo2 <- ksd2 - ksd_y
      out <- data.frame(
        ksd_y = ksd_y,
        ksd1 = ksd1,
        eo1 = eo1,
        ksd2 = ksd2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final KDD values with expected optimism penalty
  res$ksd_unpenalized[j] <- kdd(
    object1 = fit, object2 = fit,
    modified = FALSE, symmetric = TRUE
  )
  res$ksd_final1[j] <- res$ksd_unpenalized[j] + res$eo1[j]
  res$ksd_final2[j] <- res$ksd_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_ksd <- res

# Simulation for modified KDD -- KDD(R) ####
# What this means is that the model-based covariance matrix used in KDD is based
# on the working correlation structure similar to modified CIC.
set.seed(1997)
res <- data.frame(
  wcorstr = rep(x = wcorstr, each = nsims),
  ksd_unpenalized = NA,
  ksd_y = NA,
  ksd1 = NA,
  eo1 = NA,
  ksd_final1 = NA,
  ksd2 = NA,
  eo2 = NA,
  ksd_final2 = NA
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

  # eo ####
  r <- 100 # replicates

  ## Generating replicates based on the smallest model ####
  ## Smallest model is an intercept only model with an independence working
  ## correlation structure. Thus, we generate any iid N(mu, sigma), where
  ## mu and sigma are constant and fixed.
  yy <- matrix(
    data = rnorm(
      n = r * length(fit$y),
      mean = 0,
      sd = 1
    ),
    nrow = r, ncol = length(fit$y)
  )
  zz <- matrix(
    data = rnorm(
      n = r * length(fit$y),
      mean = 0,
      sd = 1
    ),
    nrow = r, ncol = length(fit$y)
  )
  parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
  ## Calculating replicate expected optimism values ####
  eo <- parallel::parLapply(
    cl = cl,
    X = seq_len(r),
    fun = function(x) {
      mod_y <- geepack::geeglm(
        formula = yy[x, ] ~ 1, id = fit$id, corstr = fit$corstr
      )
      mod_z <- geepack::geeglm(
        formula = zz[x, ] ~ 1, id = fit$id, corstr = fit$corstr
      )
      ksd_y <- kdd(
        object1 = mod_y, object2 = mod_y,
        modified = TRUE, symmetric = TRUE
      )
      ksd1 <- kdd(
        object1 = mod_y, object2 = mod_z,
        modified = TRUE, symmetric = TRUE
      )
      ksd2 <- kdd(
        object1 = mod_z, object2 = mod_y,
        modified = TRUE, symmetric = TRUE
      )
      eo1 <- ksd1 - ksd_y
      eo2 <- ksd2 - ksd_y
      out <- data.frame(
        ksd_y = ksd_y,
        ksd1 = ksd1,
        eo1 = eo1,
        ksd2 = ksd2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final KDD values with expected optimism penalty
  res$ksd_unpenalized[j] <- kdd(
    object1 = fit, object2 = fit,
    modified = TRUE, symmetric = TRUE
  )
  res$ksd_final1[j] <- res$ksd_unpenalized[j] + res$eo1[j]
  res$ksd_final2[j] <- res$ksd_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_modified_ksd <- res

# Exporting ####
save(
  res_ksd, res_modified_ksd,
  file = "./outputs/norm-ksd-eo-int-sim/sim_res.RData"
)