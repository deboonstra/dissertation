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
work_corstr <- c("independence", "exchangeable", "ar1", "unstructured")

# Simulation ####
set.seed(1997)
res_cic_corstr <- list()
res_qic_corstr <- list()
res_cic_qic <- list()
res_cic_qicu <- list()
res_qic_qic <- list()
res_qic_qicu <- list()
pb <- utils::txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr)
  y <- dat$y
  X <- dat$X
  id <- dat$id

  ## Correlation structure selection ####
  cic_corstr <- rep(NA, length(work_corstr))
  qic_corstr <- rep(NA, length(work_corstr))
  for (k in seq_along(work_corstr)) {
    fit_corstr <- geepack::geeglm(y ~ X, id = id, corstr = work_corstr[k])
    cic_corstr[k] <- cic(fit_corstr)
    qic_corstr[k] <- qic(fit_corstr)
  }

  ### Determining correlation structure ####
  ### that resulted in the minimum CIC and QIC values
  w_c_corstr <- which.min(cic_corstr)
  w_q_corstr <- which.min(qic_corstr)
  cic_best_corstr <- work_corstr[w_c_corstr]
  qic_best_corstr <- work_corstr[w_q_corstr]

  ### Updating results ####
  res_cic_corstr[[j]] <- list(
    ic = cic_corstr, work_corstr = work_corstr,
    ic_min = cic_corstr[w_c_corstr], corstr = cic_best_corstr
  )
  res_qic_corstr[[j]] <- list(
    ic = qic_corstr, work_corstr = work_corstr,
    ic_min = qic_corstr[w_q_corstr], corstr = qic_best_corstr
  )

  ## Mean structure selection ####

  ### Produce best subsets ####
  best_sub <- expand_cols(X)

  ### Fitting models ####
  c_q <- rep(NA, length(best_sub)) # cic and qic
  c_qu <- rep(NA, length(best_sub)) # cic and qicu
  q_q <- rep(NA, length(best_sub)) # qic and qic
  q_qu <- rep(NA, length(best_sub)) # qic and qicu
  penalty_c_q <- rep(NA, length(best_sub))
  penalty_c_qu <- rep(NA, length(best_sub))
  penalty_q_q <- rep(NA, length(best_sub))
  penalty_q_qu <- rep(NA, length(best_sub))
  gf_c <- rep(NA, length(best_sub)) # goodness-of-fit for cic
  gf_q <- rep(NA, length(best_sub)) # goodness-of-fit for qic
  for (i in seq_along(best_sub)) {
    XX <- X[, best_sub[[i]]]
    fit_mean_c <- geepack::geeglm(y ~ XX, id = id, corstr = cic_best_corstr)
    fit_mean_q <- geepack::geeglm(y ~ XX, id = id, corstr = qic_best_corstr)
    #### Obtaining information criterion values
    ##### CIC based
    c_q[i] <- qic(fit_mean_c)
    c_qu[i] <- qicu(fit_mean_c)
    penalty_c_q[i] <- 2 * cic(fit_mean_c)
    penalty_c_qu[i] <- 2 * nparms(fit_mean_c)
    gf_c[i] <- gof(fit_mean_c)
    ##### QIC based
    q_q[i] <- qic(fit_mean_q)
    q_qu[i] <- qicu(fit_mean_q)
    penalty_q_q[i] <- 2 * cic(fit_mean_q)
    penalty_q_qu[i] <- 2 * nparms(fit_mean_q)
    gf_q[i] <- gof(fit_mean_q)
  }
  ### Determining mean structure ####
  ### that resulted in the minimum QIC or QICu values
  w_c_q <- which.min(c_q)
  w_c_qu <- which.min(c_qu)
  w_q_q <- which.min(q_q)
  w_q_qu <- which.min(q_qu)

  ## Updating results ####
  res_cic_qic[[j]] <- list(
    ic = c_q, gof = gf_c, penalty = penalty_c_q,
    min = w_c_q, ic_min = c_q[w_c_q], vars = best_sub[[w_c_q]],
    corstr = cic_best_corstr
  )
  res_cic_qicu[[j]] <- list(
    ic = c_qu, gof = gf_c, penalty = penalty_c_qu,
    min = w_c_qu, ic_min = c_qu[w_c_qu], vars = best_sub[[w_c_qu]],
    corstr = cic_best_corstr
  )
  res_qic_qic[[j]] <- list(
    ic = q_q, gof = gf_q, penalty = penalty_q_q,
    min = w_q_q, ic_min = q_q[w_q_q], vars = best_sub[[w_q_q]],
    corstr = qic_best_corstr
  )
  res_qic_qicu[[j]] <- list(
    ic = q_qu, gof = gf_q, penalty = penalty_q_qu,
    min = w_q_qu, ic_min = q_qu[w_q_qu], vars = best_sub[[w_q_qu]],
    corstr = qic_best_corstr
  )

  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation results ####
if (!dir.exists("./outputs/both/norm-two-stage-qicu-sim/")) {
  dir.create("./outputs/both/norm-two-stage-qicu-sim/")
}

save(
  res_cic_corstr, res_qic_corstr,
  res_cic_qic, res_cic_qicu,
  res_qic_qic, res_qic_qicu,
  file = "./outputs/both/norm-two-stage-qicu-sim/norm_two_stage_qicu_sim.RData"
)