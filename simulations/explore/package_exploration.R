# The function of this script file is explore which R package will be the most
# useful as of May 2023.

# Loading libraries and functions ####
R <- list.files(path = "./simulations/R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Defining global data simulation settings ####
N <- 200
n <- 4
beta <- c(2.0, 3.0, 0.5, 0, 0, 0)
rho <- 0.5

# Simulating data ####
set.seed(1997)
dat <- sim_data(N = N, n = n, beta = beta, rho = rho)
y <- dat$y
X <- dat$X
id <- dat$id

# Model fitting ####

## geeM ####
library(Matrix)
library(geeM)
??geeM::geem
fit1 <- geem(y ~ X, id = id, corstr = "exchangeable")
summary(fit1)
fit1$alpha
fit1$biggest.R.alpha
fit1$phi

## geepack ####
library(geepack)
library(Matrix)
library(expm)
??geepack::geeglm
fit2 <- geeglm(y ~ X, id = id, corstr = "unstructured")
summary(fit2)
f2 <- fit2$geese
f2$gamma
f2$alpha
ra <- diag(1, n)
index <- lapply(
    gsub("alpha.", "", names(f2$alpha)),
    function(x) {
        hold <- unlist(strsplit(x, split = ":"))
        hold <- as.integer(hold)
        data.frame(i = hold[1], j = hold[2]) 
    }
)
index <- dplyr::bind_rows(index)
index <- dplyr::bind_rows(
    index,
    data.frame(i = index$j, j = index$i)
)
index$corr <- rep(unname(f2$alpha), 2)
for (i in 1:nrow(index)) {
    ii <- index[i, 1]
    jj <- index[i, 2]
    ra[ii, jj] <- index[i, 3]
}
V <- f2$gamma * ra
sqrt_V <- expm::sqrtm(V)
inv_sqrt_V <- Matrix::solve(sqrt_V)
inv_sqrt_V <- lapply(seq_len(N), function(x) inv_sqrt_V)
inv_sqrt_V <- Matrix::bdiag(inv_sqrt_V)
yy <- as.matrix(inv_sqrt_V %*% y)

fit2b <- geeglm(yy ~ X, id = id, corstr = "independence")
QIC(fit2b)
summary(fit2b)

## gee ####
library(gee)
??gee::gee
fit3 <- gee(y ~ X, id = id, corstr = "exchangeable")
summary(fit3)
fit3$working.correlation
fit3$scale
QIC(fit3)
