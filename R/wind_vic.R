# The function of this script file is to create a function called wind_vic that
# will identify the Occam's window interval to examine.

# N, Integer value indicating the number of clusters
# min_vic,  Numeric value indicating minimum VIC statistic under consideration
# ideal, Numeric value indicating the ideal value of criterion
# pct, Numeric value betweent 0 and 1 indicating the Gamma percentile to
# calculate. Default is 0.9

wind_vic <- function(N, min_vic, ideal, pct = 0.9) {
  # checking parameter values ####
  if (!is.numeric(N)) {
    stop("N must be an integer value.")
  }
  if (!is.numeric(min_vic)) {
    stop("min_vic must be a numeric value.")
  }
  if (!is.numeric(ideal) || ideal > min_vic) {
    stop("ideal must be a numeric value.")
  }
  if (!is.numeric(pct) || pct < 0 || pct > 1) {
    stop("pct must be a numeric value between 0 and 1.")
  }

  # Calculating window ####
  window <- qgamma(p = pct, shape = min_vic * N, rate = N) - ideal

  # Returning window
  return(window)
}