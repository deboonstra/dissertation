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
dat <- sim_data(N = N, n = n, beta = beta, rho = rho, corstr = "exchangeable")

# Transform y ####
transform_dat <- transform_y(dat)
yy <- transform_dat$y
XX <- transform_dat$X

# Produce best subsets
best_sub <- unlist(
    list(expand_cols(dat$X), expand_cols(dat$X), expand_cols(dat$X)),
    recursive = FALSE
)
wc <- rep(c("independence", "exchangeable", "ar1"), each = 63)

# Fitting independence model
qic <- rep(NA, length(best_sub))
gof <- rep(NA, length(best_sub))
penalty <- rep(NA, length(best_sub))
pb <- txtProgressBar(0, length(best_sub), style = 3)
for (i in seq_along(best_sub)) {
    xx <- dat$X[, best_sub[[i]]]
    f1 <- geepack::geeglm(yy ~ xx, id = dat$id, corstr = "independence")
    qic[i] <- unname(geepack::QIC(f1)[1])
    gof[i] <- -2 * unname(geepack::QIC(f1)[3])
    penalty[i] <- 2 * unname(geepack::QIC(f1)[4])
    setTxtProgressBar(pb, i)
    if (i == length(best_sub)) {
        close(pb)
    }
}

which.min(qic)
best_sub[[which.min(qic)]]

w <- which(beta != 0)
count <- 0
for (k in seq_len(nsims)) {
    count <- count + identical(res[[k]]$vars, w)
}