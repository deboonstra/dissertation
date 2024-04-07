# The function of this script file is to run a simulation investigating the
# relationship between the parameter estimates before and after the
# transformation of the response and design matrix.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200L
n <- 4L
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- "exchangeable"

# Simulation ####
set.seed(1997)
beta_gee <- matrix(NA, nrow = nsims, ncol = length(beta) + 1) # intercpet
beta_star <- matrix(NA, nrow = nsims, ncol = length(beta) + 1) # intercept
pb <- utils::txtProgressBar(0, nsims, style = 3)
for (j in seq_len(nsims)) {
  ## Simulating data ####
  dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = corstr)
  y <- dat$y
  X <- dat$X
  id <- dat$id
  ## Transforming data ####
  transform_data <- transform_y(data = dat, corstr = "unstructured")
  y_star <- transform_data$y
  X_star <- transform_data$X
  ## Fitting two models
  fit_gee <- geepack::geeglm(y ~ X, id = id, corstr = corstr)
  fit_star <- stats::glm(y ~ X)
  ## Obtaining coefficients
  beta_gee[j, ] <- unname(stats::coefficients(fit_gee))
  beta_star[j, ] <- unname(stats::coefficients(fit_star))
  ## Updating progress bar ####
  utils::setTxtProgressBar(pb, j)
  if (j == nsims) {
    close(pb)
  }
}

# Save simulation ####
if (!dir.exists("./outputs/meanstr/norm-glsx-beta-sim/")) {
  dir.create("./outputs/meanstr/norm-glsx-beta-sim/")
}
save(
  beta_gee, beta_star,
  file = "./outputs/meanstr/norm-glsx-beta-sim/norm_glsx_beta_sim.RData"
)