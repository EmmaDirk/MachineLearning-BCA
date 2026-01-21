# This script contains a function to sample the coefficients with which time-invariant confounders (e.g. age),
# affect the variables X and Y at time t = 1. 
# Since the autoregessive effects are denoted using beta, and the cross-lagged effects are denoted using gamma,
# we denote the confounder effects using delta. As such, here we sample the delta coefficients at time t = 1. 
#
# The function samples the coefficients such that:
# - the total variance explained by the confounders at time t = 1 is R2_total
# - the fraction of variance explained by the non-linear confounders is R2_nonlin
#
# If the non-linear confounders explain any variance (>0), we include all interaction terms of order 2..k
# where k is the number of confounders. E.g. for k = 3 confounders c1, c2, c3, we include the terms:
# - c1:c2
# - c1:c3
# - c2:c3
# - c1:c2:c3
#
# This is accomplished by sampling random coefficients from a normal distribution, and scaling them
# such that the sum of squares of the coefficients equals the desired R^2 values for linear and non-linear parts.
# The function checks whether the sum of squares of the sampled coefficients is within the desired range. 
# If not, it resamples until valid coefficients are found, or the maximum number of tries is reached.
# -----------------------------------------------------------------------------------------------------------------

sample_delta_1 <- function(
  k,                                                             # number of confounders
  R2_total,                                                      # total variance explained by the confounders at t=1
  R2_nonlin = 0.0,                                               # fraction of R2_total allocated to non-linear terms
  min_abs   = 0.01,                                              # minimum absolute value for each delta coefficient
  max_abs   = 0.60,                                              # maximum absolute value for each delta coefficient
  max_tries = 100000                                             # maximum sampling attempts
) {

  # check R2_nonlin is a valid fraction
  if (R2_nonlin < 0 || R2_nonlin > 1)
    stop("R2_nonlin must be between 0 and 1")

  # number of linear terms (always k main effects)
  k_lin <- k

  # number of non-linear terms (all interaction terms of order 2..k)
  # if R2_nonlin > 0
  k_non <- if (R2_nonlin > 0) {

    # calculate number of non-linear terms as the sum of combinations of k taken m at a time for m = 2..k
    sum(sapply(2:k, function(m) choose(k, m)))

    # otherwise, the number of non-linear terms is 0
  } else 0

  # split total R^2 into linear and non-linear parts
  # of the total variance explained by confounders (R2_total), 1 - R2_nonlin is allocated to linear terms
  R2_lin <- (1 - R2_nonlin) * R2_total

  # of the total variance explained by confounders (R2_total), R2_nonlin is allocated to non-linear terms
  R2_non <- R2_nonlin * R2_total

  # start the loop: for i in the max number of tries
  # max tries is set because we may not be able to sample valid coefficients that meet the min/max abs criteria
  for (i in seq_len(max_tries)) {

    # sample the delta coefficients for the linear effects of the confounders on X
    # if we have linear terms (i.e. k > 0) and R2_lin > 0, then:
    dx_lin <- if (k_lin > 0 && R2_lin > 0) {

      # sample k random numbers from a normal distribution
      v <- rnorm(k_lin)

      # normalize so the sum of squares of this vector is 1,
      # then scale so the sum of squares is R2_lin
      sqrt(R2_lin) * v / sqrt(sum(v^2))

      # otherwise, no linear effects
    } else numeric(0)  

    # sample the delta coefficients for the linear effects of the confounders on Y
    # if we have linear terms (i.e. k > 0) and R2_lin > 0, then:
    dy_lin <- if (k_lin > 0 && R2_lin > 0) {

      # sample k random numbers from a normal distribution
      v <- rnorm(k_lin)

      # normalize + scale to hit R2_lin
      sqrt(R2_lin) * v / sqrt(sum(v^2))

      # otherwise, no linear effects
    } else numeric(0)

    # initialize non-linear vectors as empty
    dx_non <- dy_non <- numeric(0)

    # sample the delta coefficients for the non-linear effects of the confounders on X and Y
    # if we have non-linear terms (i.e. k_non > 0) and R2_non > 0, then:
    if (k_non > 0 && R2_non > 0) {

      # sample k_non random numbers from a normal distribution for X
      v <- rnorm(k_non)

      # normalize + scale to hit R2_non
      dx_non <- sqrt(R2_non) * v / sqrt(sum(v^2))

      # sample k_non random numbers from a normal distribution for Y
      v <- rnorm(k_non)

      # normalize + scale to hit R2_non
      dy_non <- sqrt(R2_non) * v / sqrt(sum(v^2))
    }

    # combine linear + non-linear parts into full delta vectors
    dx <- c(dx_lin, dx_non)
    dy <- c(dy_lin, dy_non)

    # check if all absolute values are within the specified bounds
    if (all(abs(c(dx, dy)) >= min_abs & abs(c(dx, dy)) <= max_abs)) {

      # if so, return the Delta matrix (2 rows: X and Y)
      D1 <- rbind(dx, dy)

      # set row names
      rownames(D1) <- c("X", "Y")

      # set column names:
      # - first k are main effects: c1..ck
      # - then (if present) interaction names like c1:2, c1:3, ..., c1:2:3, ...
      if (k_non == 0) {

        # only linear terms
        colnames(D1) <- paste0("c", 1:k)

      } else {

        # build interaction names in the same order as combn(1:k, m)
        int_names <- unlist(
          lapply(2:k, function(m)
            combn(1:k, m, FUN = function(ix)
              paste0("c", paste(ix, collapse=":")))
          )
        )

        # concatenate main effect names + interaction names
        colnames(D1) <- c(paste0("c", 1:k), int_names)
      }

      return(D1)
    }
  }

  # if the loop doesn't return a valid Delta matrix, throw an error
  stop("Failed to sample Delta at t=1 within max_tries")
}
