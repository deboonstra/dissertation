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
# strict_upper, Logical indicating if a strict upper bound should be used. If
# TRUE then only the values below the upper bound are considered and the
# resulting implementation of Occam's window could select none of the
# structures. If FALSE, then the structure with minimum value will be selected
# if all of the values fall outside of the Occam's window.

# Details: If length(window) ==  2 then the value of x closet in absolute value
# to ideal is choosen.

occam <- function(
    x, ideal, window = 1, nparms = NULL, mv = 1, strict_upper = TRUE) {
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
  if (!is.logical(strict_upper)) {
    stop("strict_upper must be a logical.")
  }

  # Finding the number of parameters

  ## mv ####
  if (sum(names(x) %in% mv_names) > 0) {
    mv <- mv
  }

  ##  Number of parameters estimated ####
  if (!is.null(nparms)) {
    nparms <- nparms
  } else {
    nparms <- rep(NA, times = length(x))
    ### Independence ####
    if ("independence" %in% names(x)) {
      w_nparms <- which(names(x) == "independence")
      nparms[w_nparms] <- 0
    }

    ### Exchangeable ####
    if ("exchangeable" %in% names(x)) {
      w_nparms <- which(names(x) == "exchangeable")
      nparms[w_nparms] <- 1
    }

    ### AR-M, Stationary and non-statationary M-dependent ####
    if (sum(names(x) %in% mv_names) > 0) {
      w_nparms <- which(names(x) %in% mv_names)
      nparms[w_nparms] <- mv[w_nparms]
    }

    ### Unstructured ####
    if ("unstructured" %in% names(x)) {
      w_nparms <- which(names(x) == "unstructured")
      nparms[w_nparms] <- (ideal * (ideal - 1)) / 2
    }
  }


  # Implementing Occam's window ####

  ## Determining which values within window ####
  if (length(window) == 1) {
    w <- which(x <= (ideal + window))
  } else if (length(window) == 2) {
    if (window[1] > 0) {
      window[1] <- -1 * window[1]
    }
    if (window[2] < 0) {
      window[2] <- -1 * window[2]
    }
    w <- which(x >= ideal + window[1] & x <= ideal + window[2])
  }

  ## Selecting value with minimum number of parameters ####
  if (length(w) > 0) {
    ### Selecting values within window ####
    x_wind <- x[w]
    nparms_wind <- nparms[w]

    ### Selecting the value with the minimum number of parameters ####
    w_sel <- which(nparms_wind == min(nparms_wind))
    if (length(w_sel) == 1) {
      x_sel <- x_wind[w_sel]
      nparms_sel <- nparms_wind[w_sel]
    } else if (length(w_sel) > 1) {
      if (length(window) == 1) {
        x_sel <- x_wind[w_sel]
        nparms_sel <- nparms_wind[w_sel]
        x_sel <- x_sel[which(x_sel == min(x_sel))]
        nparms_sel <- nparms_sel[which(x_sel == min(x_sel))]
      } else if (length(window) == 2) {
        x_sel <- x_wind[w_sel]
        nparms_sel <- nparms_wind[w_sel]
        abs_diff <- abs(x_sel - ideal)
        x_sel <- x_sel[which(abs_diff == min(abs_diff))]
        nparms_sel <- nparms_sel[which(abs_diff == min(abs_diff))]
      }
    }

    ### Creating output data frame given at least one value in window ####
    out <- data.frame(
      name = names(x_sel),
      value = unname(x_sel),
      ideal = unname(ideal),
      nparms = unname(nparms_sel)
    )
  } else {
    ### Creating output data frame given no value in window and
    ### strict upper is TRUE resulting in NO structure selected ####
    if (strict_upper) {
      out <- data.frame(
        name = NA,
        value = NA,
        ideal = unname(ideal),
        nparms = NA
      )
    } else {
      ### Creating output data frame given no value in window and ####
      ### strict upper bound is FALSE resulting in minimum value selected
      out <- data.frame(
        name = names(x)[which(x == min(x))],
        value = x[which(x == min(x))],
        ideal = unname(ideal),
        nparms = nparms[which(x == min(x))]
      )
    }
  }

  # Returning data frame ####
  return(out)
}