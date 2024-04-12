# The function of this script file is to explore the null distrubtion of the
# delta statistic. The null distribution will be characterized by the delta
# values produced by the working correlation structures that are also the
# generating correlation structures.

# Loading libraries and functions ####
R <- list.files(path = "./R", pattern = "*.R", full.names = TRUE)
sapply(R, source, .GlobalEnv)

# Creating output sub-directory ####
if (!dir.exists("./outputs/corstr/norm-delta-sim/norm-delta-null-dist/")) {
  dir.create("./outputs/corstr/norm-delta-sim/norm-delta-null-dist/")
}

# Defining simulation information ####
N <- c(rep(200, 7), 300, 240, 150, 120, 100, 80, 75)
n <- c(4, 5, 6, 7, 8, 9, 10, 4, 5, 8, 10, 12, 15, 16)
corstr <- c("exchangeable", "ar1")
dim_beta <- c(3, 5, 3, 5, 3, 5)
alpha <- c(0.50, 0.50, 0.05, 0.05, 0.70, 0.70)

# Importing simulation results ####

## List of simulation results ####
res_list <- list.files(
  path = "./outputs/corstr/norm-delta-sim",
  pattern = ".rds",
  recursive = TRUE,
  full.names = TRUE
)

### Removing null distribution results if this is ran again ####
exclude <- c(
  paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-null-dist/norm_delta_null_dist.rds"
  ),
  paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-null-dist/norm_delta_realized_values.rds"
  )
)
if (all(file.exists(exclude))) {
  res_list <- res_list[!(res_list %in% exclude)]
}

## Importing ####
res <- lapply(X = res_list, FUN = readRDS)
names(res) <- substr(
  x = res_list,
  start = regexpr(pattern = "delta_", text = res_list, perl = TRUE),
  stop = regexpr(pattern = "_sim.rds", text = res_list, perl = TRUE) - 1
)

# Finding the null distribution ####

## Creating a basis data.frame to store all the results ####
null_dist <- data.frame(
  simulation = rep(x = names(res), each = length(N) * length(corstr)),
  dim_beta = rep(x = dim_beta, each = length(N) * length(corstr)),
  alpha = rep(x = alpha, each = length(N) * length(corstr)),
  N = rep(x = N, times = length(res) * length(corstr)),
  n = rep(x = n, times = length(res) * length(corstr)),
  corstr = rep(x = rep(x = corstr, each = length(N)), times = length(res)),
  shape = rep(x = NA, times = length(res) * length(N) * length(corstr)),
  rate = rep(x = NA, times = length(res) * length(N) * length(corstr)),
  scale = rep(x = NA, times = length(res) * length(N) * length(corstr)),
  mean = rep(x = NA, times = length(res) * length(N) * length(corstr)),
  var = rep(x = NA, times = length(res) * length(N) * length(corstr))
)

## Organizing results to match basis data.frame ####
for (k in seq_along(res)) {
  res[[k]]$simulation <- names(res)[k]
  res_cs <- subset(
    x = res[[k]],
    subset = corstr == "exchangeable",
    select = c(simulation, N, n, corstr, exchangeable)
  )
  res_ar1 <- subset(
    x = res[[k]],
    subset = corstr == "ar1",
    select = c(simulation, N, n, corstr, ar1)
  )
  colnames(res_cs) <- c("simulation", "N", "n", "corstr", "delta")
  colnames(res_ar1) <- c("simulation", "N", "n", "corstr", "delta")
  res[[k]] <- rbind.data.frame(res_cs, res_ar1)
  res[[k]]$dim_beta <- rep(x = dim_beta[k], times = nrow(res[[k]]))
  res[[k]]$alpha <- rep(x = alpha[k],  times = nrow(res[[k]]))
  res[[k]] <- subset(
    x = res[[k]],
    select = c(simulation, dim_beta, alpha, N, n, corstr, delta)
  )
}
res <- dplyr::bind_rows(res)

## Find the parameters of the null distribution ####
for (j in seq_len(nrow(null_dist))) {
  # Pulling the delta values related to the basis data.frame ####
  x <- subset(
    x = res,
    subset = (
      simulation == null_dist$simulation[j] &
        N == null_dist$N[j] &
        n == null_dist$n[j] &
        corstr == null_dist$corstr[j]
    )
  )
  x <- c(x$delta)

  # Find the distribution parameters ####
  dist_info <- MASS::fitdistr(x = x, densfun = "gamma")

  # Filling the basis data.frame ####
  null_dist$shape[j] <- dist_info$estimate[1]
  null_dist$rate[j] <- dist_info$estimate[2]
  null_dist$scale[j] <- 1 / null_dist$rate[j]
  null_dist$mean[j] <- null_dist$shape[j] / null_dist$rate[j]
  null_dist$var[j] <- null_dist$shape[j] / (null_dist$rate[j]^2)
}

# Exporting results ####
saveRDS(
  object = null_dist,
  file = paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-null-dist/norm_delta_null_dist.rds"
  )
)
saveRDS(
  object = res,
  file = paste0(
    "./outputs/corstr/norm-delta-sim/",
    "norm-delta-null-dist/norm_delta_realized_values.rds"
  )
)