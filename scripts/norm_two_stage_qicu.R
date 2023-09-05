# The function of this script file is to run a simulation to investigate the
# selection properties of a two-stage approach, where the correlation structure
# is determined first and used to select the mean structure.
# See 2023_08_31 notes for more information.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- "exchangeable"
work_corstr <- c("exchangeable", "ar1", "unstructured")

# Simulation ####
set.seed(1997)
res_qic <- list()
res_qicu <- list()
res_cic <- list()
pb <- utils::txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = n, n = n, beta = beta, rho = rho, corstr = corstr)
  y <- dat$y
  X <- dat$X
  id <- dat$id

  ## Correlation structure selection ####
  cic_corstr <- rep(NA, length(work_corstr))
  for (k in seq_along(work_corstr)) {
    fit_corstr <- geepack::geeglm(y ~ X, id = id, corstr = work_corstr[k])
    cic_corstr[k] <- cic(fit_corstr)
  }

  ### Determining correlation structure ####
  ### that resulted in the minimum CIC value
  w_c <- which.min(cic_corstr)
  best_corstr <- work_corstr[w_c]

  ### Updating results ####
  res_cic[[j]] <- list(
    ic = cic_corstr, ic_min = cic_corstr[w_c], corstr = best_corstr
  )

  ## Mean structure selection ####

  ### Produce best subsets ####
  best_sub <- expand_cols(X)

  ### Fitting models ####
  q <- rep(NA, length(best_sub)) # qic
  qu <- rep(NA, length(best_sub)) # qicu
  penalty_q <- rep(NA, length(best_sub))
  penalty_qu <- rep(NA, length(best_sub))
  gf <- rep(NA, length(best_sub)) # goodness-of-fit
  for (i in seq_along(best_sub)) {
    XX <- X[, best_sub[[i]]]
    fit_mean <- geepack::geeglm(y ~ XX, id = id, corstr = best_corstr)
    q[i] <- qic(fit_mean)
    qu[i] <- qicu(fit_mean)
    penalty_q[i] <- 2 * cic(fit_mean)
    penalty_qu[i] <- 2 * nparms(fit_mean)
    gf[i] <- gof(fit_mean)
  }
  ### Determining mean structure ####
  ### that resulted in the minimum QIC or QICu values
  w_q <- which.min(q)
  w_qu <- which.min(qu)

  ## Updating results ####
  res_qic[[j]] <- list(
    ic = qic, gof = gf, penalty = penalty_q,
    min = w_q, ic_min = q[w_q], vars = best_sub[[w_q]]
  )
  res_qicu[[j]] <- list(
    ic = qicu, gof = gf, penalty = penalty_qu,
    min = w_qu, ic_min = qu[w_qu], vars = best_sub[[w_qu]]
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation results ####
if (!dir.exists("./outputs/norm_two_stage_qicu")) {
  dir.create("./outputs/norm_two_stage_qicu")
}

save(
  res_cic, res_qic, res_qicu,
  file = "./outputs/norm_two_stage_qicu/norm_two_stage_qicu.RData"
)