# The function of this script file is to run a simulation study where the
# performance of AIC and QIC will be compared. This comparison can only occur
# if the working correlation structure is independence.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0
corstr <- "independence"

# Simulation ####
set.seed(1997)
res_qic <- list()
res_aic <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr)
  y <- dat$y
  X <- dat$X
  id <- dat$id
  ## Produce best subsets ####
  best_sub <- expand_cols(X)
  # Fitting models ####
  q <- gf_q <- pen_q <- rep(NA, length(best_sub))
  a <- gf_a <- pen_a <- rep(NA, length(best_sub))
  for (i in seq_along(best_sub)) {
    XX <- X[, best_sub[[i]]]
    f_q <- geepack::geeglm(y ~ XX, id = id, corstr = corstr)
    f_a <- stats::glm(y ~ XX)
    q[i] <- qic(f_q)
    gf_q[i] <- gof(f_q)
    pen_q[i] <- 2 * cic(f_q)
    a[i] <- stats::AIC(f_a)
    gf_a[i] <- -2 * logLik(f_a)
    pen_a[i] <- a[i] - gf_a[i]
  }
  # Determining the best_sub index ####
  # that resulted in the minimum AIC/QIC values
  w_q <- which.min(q)
  w_a <- which.min(a)
  # Updating res objects ####
  res_qic[[j]] <- list(
    ic = q, gof = gf_q, penalty = pen_q,
    min = w_q, ic_min = q[w_q], vars = best_sub[[w_q]]
  )
  res_aic[[j]] <- list(
    ic = a, gof = gf_a, penalty = pen_a,
    min = w_a, ic_min = a[w_a], vars = best_sub[[w_a]]
  )
  setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulations ####
save(
  res_qic, res_aic,
  file = "./outputs/norm-aic-sim/norm_aic_sim.RData"
)