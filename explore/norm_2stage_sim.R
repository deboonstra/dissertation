# The function of this script file to run a simulation study exploring using CIC
# and QIC in a two-stage approach to select the proper mean and correlation
# structures.

# Loading libraries and functions ####
R <- list.files(path = "./simulations/R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
nsims <- 100L
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5
corstr <- c("independence", "exchangeable", "ar1", "userdefined")

# Simulation ####
set.seed(1997)
corstr_true <- rep(corstr, each = nsims)
mat <- un(n = n)
res <- rep(NA, length(corstr_true))
pb <- txtProgressBar(0, length(corstr_true), style = 3)
for (j in seq_along(corstr_true)) {
  ## Simulating data ####
  dat <- sim_data(
    N = N, n = n, beta = beta, rho = rho, corstr = corstr_true[j], zcor = mat
  )
  y <- dat$y
  X <- dat$X
  id <- dat$id
  ## Selecting correlation structure ####
  ### Obtaining CIC values ####
  ### based on saturated model
  corstr_ic <- sapply(
    X = corstr,
    FUN = function(x) {
      if (x != "userdefined") {
        cic(geepack::geeglm(y ~ X, id = id, corstr = x))
      } else {
        cic(geepack::geeglm(y ~ X, id = id, corstr = "unstructured"))
      }
    }
  )
  ### Selecting corstr with mininum CIC ####
  corstr_w <- names(corstr_ic)[which.min(corstr_ic)]
  corstr_w <- ifelse(corstr_w == "userdefined", "unstructured", corstr_w)
  res[j] <- corstr_w
  setTxtProgressBar(pb, j)
  if (j == length(corstr_true)) {
    close(pb)
  }
}