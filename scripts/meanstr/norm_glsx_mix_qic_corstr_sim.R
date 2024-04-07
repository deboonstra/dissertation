# The function of this script file is to run a simulation to determine the
# effectiveness of using a mixed GLS approach with GEEs. The GOF term will be
# a by-product of the transformed y and X; however, the penalty term will be
# based on the original data.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5

# Simulation ####
set.seed(1997)
res_quasi <- list()
res_lik <- list()
pb <- txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
  y <- dat$y
  X <- dat$X
  id <- dat$id
  ## Transform data ####
  transform_data <- transform_y(dat)
  yt <- transform_data$y
  Xt <- transform_data$X
  ## Produce best subsets ####
  bs_temp <- expand_cols(dat$X)
  best_sub <- unlist(
    list(bs_temp, bs_temp, bs_temp, bs_temp),
    recursive = FALSE
  )
  wc <- rep(
    x = c("independence", "exchangeable", "ar1", "unstructured"),
    each = length(bs_temp)
  )
  ## Fitting models ####
  ic_quasi <- ic_lik <- gf_quasi <- gf_lik <- pen <- corstr <- rep(NA, length(best_sub))
  for (i in seq_along(best_sub)) {
    ### Pulling covariate matrices ####
    XX <- X[, best_sub[[i]]]
    XXt <- Xt[, best_sub[[i]]]
    ### Fitting original data for penalty ####
    f <- geepack::geeglm((y ~ XX), id = id, corstr = wc[i])
    ### Fitting transformed data for GOF ####
    #### Quasi-likelihood
    ft_quasi <- geepack::geeglm(yt ~ XXt, id = id, corstr = "independence")
    #### Likelihood
    ft_lik <- stats::glm(yt ~ XXt)
    ### Obtaining information criterion values ####
    gf_quasi[i] <- gof(ft_quasi)
    gf_lik[i] <- -2 * stats::logLik(ft_lik)
    pen[i] <- 2 * cic(f)
    ic_quasi[i] <- gf_quasi[i] + pen[i]
    ic_lik[i] <- gf_lik[i] + pen[i]
  }
  ## Determining the best_sub index ####
  ## that resulted in the minimum IC
  w_quasi <- which.min(ic_quasi)
  w_lik <- which.min(ic_lik)
  ## Updating res object(s) ####
  res_quasi[[j]] <- list(
    ic = ic_quasi, gof = gf_quasi, penalty = pen, corstr = wc,
    min = w_quasi, ic_min = ic_quasi[w_quasi],
    corstr_min = wc[w_quasi], vars = best_sub[[w_quasi]]
  )
  res_lik[[j]] <- list(
    ic = ic_lik, gof = gf_lik, penalty = pen, corstr = wc,
    min = w_lik, ic_min = ic_lik[w_lik],
    corstr_min = wc[w_lik], vars = best_sub[[w_lik]]
  )
  setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation ####
if (!dir.exists("./outputs/norm-glsx-mix-qic-corstr-sim/")) {
  dir.create("./outputs/norm-glsx-mix-qic-corstr-sim/")
}
save(
  res_quasi, res_lik,
  file = "./outputs/norm-glsx-mix-qic-corstr-sim/norm_glsx_mix_qic_corstr_sim.RData"
)