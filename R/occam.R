# The function of this script file is to create a function that implements
# Occam's window on a name vector of values.

# This function is designed for correlation structure selection based on the
# gee package.

# x, Named numeric vector of values to select from.
# ideal, Numeric value indicating the ideal value of criterion
# window, How big the view point around the minimum should be. If length(window)
# is 1 then search is ideal + window. If the length(window) is 2 then
# the search is ideal + window[1] and ideal + window[2], where window[1] should
# be a negative value.
# nparms, numeric vector of the number of parameters estimated for each value in
# x. Thus, nparms must be same length as x if specified. Default is NULL, as
# this may be obtained from the names of x if x is named based on the corstr
# options pre
# mv, When x contains "stat_M_dep", "non_stat_M_dep", or "AR-M" then mv must be
# specified and the same length as x.

# Details: If length(window) ==  2 then the value of x closet in absolute value
# to ideal is choosen.

occam <- function(x, ideal, window = 1, nparms = NULL, mv = 1) {
  # Checking parameter values ####
  if (!is.numeric(x) || is.null(names(x))) {
    stop("x must be a named numeric vector.")
  }
  if (!is.numeric(ideal) || length(ideal) != 1) {
    stop("ideal must be a single numeric value.")
  }
  if (!is.numeric(window)) {
    stop("window must be a numeric vector.")
  }
  if (length(window) < 0 || length(window) > 2) {
    stop("window must be a numeric vector of length 1 or 2.")
  }
  if (!is.null(nparms)) {
    if (!is.numeric(nparms) || length(nparms) != length(x)) {
      stop("nparms must be a numeric vector the same length as x.")
    }
  }
  mv_names <- c("stat_M_dep", "non_stat_M_dep", "AR-M")
  if (sum(names(x) %in% mv_names) > 0) {
    if (length(mv) != length(x)) {
      stop("mv must be a numeric vector the same length as x.")
    }
  }

  # Implementing Occam's window ####

  ## Subsetting x for the values within window ####

  ### Determining which values within window ####
  if (length(window) == 1) {
    w <- which(x >= ideal & x <= (ideal + window))
  } else if (length(window) == 2) {
    if (window[1] > 0) {
      window[1] <- -1 * window[1]
    }
    if (window[2] < 0) {
      window[2] <- -1 * window[2]
    }
    w <- which(x >= ideal + window[1] & x <= ideal + window[2])
  }

  ### Subsetting ####

  #### Named vector of values to select from ####
  x_wind <- x[w]

  #### Number of parameters estimated ####
  if (!is.null(nparms)) {
    nparms <- nparms[w]
  }

  #### mv ####
  if (sum(names(x) %in% mv_names) > 0) {
    mv <- mv[w]
  }


  ## Calculating number of parameters estimated ####
  if (is.null(nparms)) {
    ### Independence ####
    if ("independence" %in% names(x_wind)) {
      w_nparms <- which(names(x_wind) == "independence")
      nparms[w_nparms] <- 0
    }

    ### Exchangeable ####
    if ("exchangeable" %in% names(x_wind)) {
      w_nparms <- which(names(x_wind) == "exchangeable")
      nparms[w_nparms] <- 1
    }

    ### AR-M, Stationary and non-statationary M-dependent ####
    if (sum(names(x_wind) %in% mv_names) > 0) {
      w_nparms <- which(names(x_wind) %in% mv_names)
      nparms[w_nparms] <- mv[w_nparms]
    }

    ### Unstructured ####
    if ("unstructured" %in% names(x_wind)) {
      w_nparms <- which(names(x_wind) == "unstructured")
      nparms[w_nparms] <- (ideal * (ideal - 1)) / 2
    }
  }

  ## Selecting value with minimum number of parameters ####
  w_sel <- which(nparms == min(nparms))
  if (length(w_sel) == 1) {
    x_sel <- x_wind[w_sel]
    nparms_sel <- nparms[w_sel]
  } else if (length(w_sel) > 1) {
    if (length(window) == 1) {
      x_sel <- x_wind[w_sel]
      nparms_sel <- nparms[w_sel]
      x_sel <- x_sel[which(x_sel == min(x_sel))]
      nparms_sel <- nparms_sel[which(x_sel == min(x_sel))]
    } else if (length(window) == 2) {
      x_sel <- x_wind[w_sel]
      nparms_sel <- nparms[w_sel]
      abs_diff <- abs(x_sel - ideal)
      x_sel <- x_sel[which.min(abs_diff)]
      nparms_sel <- nparms_sel[which.min(abs_diff)]
    }
  }

  # Creating output data frame ####
  out <- data.frame(
    name = names(x_sel),
    value = unname(x_sel),
    ideal = unname(ideal),
    nparms = unname(nparms_sel)
  )

  # Returning data frame ####
  return(out)
}