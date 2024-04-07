# The function of this script file is to run a simulation study, where the
# performance of QIC, QICu, and BQICu will be compared through the traditional
# application of GEEs and information criterions (i.e., no transformations).
# See 2023_05_2023 notes for more detail.
# These simulations were updated to include selection of correlation structures
# and including unstructured correlation structure as an option.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- c("independence", "exchangeable", "ar1", "unstructured")

# Simulation ####
set.seed(1997)
res_qic <- list()
res_qicu <- list()
res_bqicu <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
  y <- dat$y
  X <- dat$X
  id <- dat$id
  ## Produce best subsets
  bs_temp <- expand_cols(X)
  best_sub <- unlist(
    list(bs_temp, bs_temp, bs_temp, bs_temp),
    recursive = FALSE
  )
  wc <- rep(corstr, each = length(bs_temp))
  ## Fitting model ####
  q <- gf_q <- pen_q <- rep(NA, length(best_sub))
  qu <- gf_qu <- pen_qu <- rep(NA, length(best_sub))
  bqu <- gf_bqu <- pen_bqu <- rep(NA, length(best_sub))
  for (i in seq_along(best_sub)) {
    ### Pulling covariate matrices ####
    XX <- X[, best_sub[[i]]]
    ### Fitting each best subset model ####
    f <- geepack::geeglm(y ~ XX, id = id, corstr = wc[i])
    ### Obtaining information criterion values ####
    q[i] <- qic(f)
    gf_q[i] <- gof(f)
    pen_q[i] <- 2 * cic(f)
    qu[i] <- qicu(f)
    gf_qu[i] <- gof(f)
    pen_qu[i] <- qu[i] - gf_qu[i]
    bqu[i] <- bqicu(f)
    gf_bqu[i] <- gof(f)
    pen_bqu[i] <- bqu[i] - gf_bqu[i]
  }
  ## Determining the best_sub index ####
  ## that resulted in the minimum QIC,
  ## QICu, and BQICu values
  w_q <- which.min(q)
  w_qu <- which.min(qu)
  w_bqu <- which.min(bqu)
  ## Updating res object(s) ####
  res_qic[[j]] <- list(
    ic = q, gof = gf_q, penalty = pen_q, corstr = wc,
    min = w_q, ic_min = q[w_q], corstr_min = wc[w_q],
    vars = best_sub[[w_q]]
  )
  res_qicu[[j]] <- list(
    ic = qu, gof = gf_qu, penalty = pen_qu, corstr = wc,
    min = w_qu, ic_min = q[w_qu], corstr_min = wc[w_qu],
    vars = best_sub[[w_qu]]
  )
  res_bqicu[[j]] <- list(
    ic = bqu, gof = gf_bqu, penalty = pen_bqu, corstr = wc,
    min = w_bqu, ic_min = q[w_bqu], corstr_min = wc[w_bqu],
    vars = best_sub[[w_bqu]]
  )
  ## Updating progress bar ####
  setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulations ####
save(
  res_qic, res_qicu, res_bqicu,
  file = "./outputs/both/norm-bqicu-sim/norm_bqicu_sim.RData"
)