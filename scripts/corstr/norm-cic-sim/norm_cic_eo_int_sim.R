# The function of this script is to run a simulation investigating the selection
# properties of CIC penalized by expected optimism (eo). Additionally, we will
# be determining how to properly specify eo in these simulations. The mean
# difference component will not be calculated in these simulations as the
# limiting distributions of all the models will have the same mean vector.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output directory if one does not exists
if (!dir.exists("./outputs/corstr/norm-cic-sim/norm-cic-eo-int-sim/")) {
  dir.create("./outputs/corstr/norm-cic-sim/norm-cic-eo-int-sim/")
}

# Defining global data simulation settings ####
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
nsims <- 1000
wcorstr <- c("exchangeable", "ar1", "independence", "unstructured")

# Simulation for original CIC -- CIC(I) ####
# What this means is that the model-based covariance matrix used in CIC is based
# on the independence working correlation structure similar to QIC and CIC.
set.seed(1997)
res <- data.frame(
  wcorstr = rep(x = wcorstr, each = nsims),
  cic_unpenalized = NA,
  cic_y = NA,
  cic1 = NA,
  eo1 = NA,
  cic_final1 = NA,
  cic2 = NA,
  eo2 = NA,
  cic_final2 = NA
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
      cic_y <- cic(object1 = mod_y, object2 = mod_y)
      cic1 <- cic(object1 = mod_y, object2 = mod_z)
      cic2 <- cic(object1 = mod_z, object2 = mod_y)
      eo1 <- cic1 - cic_y
      eo2 <- cic2 - cic_y
      out <- data.frame(
        cic_y = cic_y,
        cic1 = cic1,
        eo1 = eo1,
        cic2 = cic2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final cic values with expected optimism penalty
  res$cic_unpenalized[j] <- cic(object1 = fit, object2 = fit)
  res$cic_final1[j] <- res$cic_unpenalized[j] + res$eo1[j]
  res$cic_final2[j] <- res$cic_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_cic <- res

# Simulation for modified CIC -- CIC(R) ####
# What this means is that the model-based covariance matrix used in CIC is based
# on the working correlation structure similar to modified CIC.
set.seed(1997)
res <- data.frame(
  wcorstr = rep(x = wcorstr, each = nsims),
  cic_unpenalized = NA,
  cic_y = NA,
  cic1 = NA,
  eo1 = NA,
  cic_final1 = NA,
  cic2 = NA,
  eo2 = NA,
  cic_final2 = NA
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
      cic_y <- cic(object1 = mod_y, object2 = mod_y, modified = TRUE)
      cic1 <- cic(object1 = mod_y, object2 = mod_z, modified = TRUE)
      cic2 <- cic(object1 = mod_z, object2 = mod_y, modified = TRUE)
      eo1 <- cic1 - cic_y
      eo2 <- cic2 - cic_y
      out <- data.frame(
        cic_y = cic_y,
        cic1 = cic1,
        eo1 = eo1,
        cic2 = cic2,
        eo2 = eo2
      )
      return(out)
    }
  )
  ## Calculating expected optimism for each simulation
  res[j, c(3, 4, 5, 7, 8)] <- colMeans(dplyr::bind_rows(eo))
  ## Calculating final cic values with expected optimism penalty
  res$cic_unpenalized[j] <- cic(object1 = fit, object2 = fit, modified = TRUE)
  res$cic_final1[j] <- res$cic_unpenalized[j] + res$eo1[j]
  res$cic_final2[j] <- res$cic_unpenalized[j] + res$eo2[j]
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nrow(res)) {
    close(pb)
    parallel::stopCluster(cl)
  }
}
res_modified_cic <- res

# Exporting ####
save(
  res_cic, res_modified_cic,
  file = "./outputs/corstr/norm-cic-sim/norm-cic-eo-int-sim/sim_res.RData"
)