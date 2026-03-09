# This script contains a function to sample the coefficients with which time-invariant confounders (e.g. age),
# affect the variables X and Y at time t = 1.
# Since the autoregessive effects are denoted using beta, and the cross-lagged effects are denoted using gamma,
# we denote the confounder effects using delta. As such, here we sample the delta coefficients at time t = 1.
#
# The function samples the coefficients such that:
# - the total variance explained by the confounders at time t = 1 is R2_total
# - the fraction of variance explained by interaction terms is R2_interaction
#
# The argument k refers to the number of base confounders (main effects): c1..ck.
# Interaction terms are only included when requested via toggles:
# - include_2way: include every possible two-way interaction (e.g. c1:c2, c1:c3, ...)
# - include_3way: include every possible three-way interaction (e.g. c1:c2:c3, c1:c2:c4, ...)
#
# The interaction variance is allocated across whichever interaction terms are included by these toggles.
# If interaction terms are requested, but the interaction variance is zero, the function throws an error.
# Likewise, if interaction variance is positive but no interaction terms are requested, the function throws an error.
#
# This is accomplished by sampling random coefficients from a normal distribution, and scaling them
# such that the sum of squares of the coefficients equals the desired R^2 values for linear and interaction parts.
# The function checks whether the absolute values of the sampled coefficients are within the desired range.
# If not, it resamples until valid coefficients are found, or the maximum number of tries is reached.
#
# This requires the following assumptions / restictions:
# - All main effects and X and Y have an expected value of 0
# - All main effects and X and Y have a variance of 1
# - All 'base' confounders are independent
# -----------------------------------------------------------------------------------------------------------------

sample_delta_1 <- function(
  k,                                                             # number of base confounders
  R2_total,                                                      # total variance explained by the confounders at t=1
  R2_interaction = 0,                                            # fraction of R2_total allocated to interaction terms
  include_2way = FALSE,                                          # include every possible two-way interaction
  include_3way = FALSE,                                          # include every possible three-way interaction
  min_abs   = 0,                                                 # minimum absolute value for each delta coefficient
  max_abs   = 1,                                                 # maximum absolute value for each delta coefficient
  max_tries = 100000                                             # maximum sampling attempts
) {

  # check R2_total is a valid fraction
  if (R2_total < 0 || R2_total > 1)
    stop("R2_total must be between 0 and 1")
  
  # check R2_interaction is a valid fraction
  if (R2_interaction < 0 || R2_interaction > 1)
    stop("R2_interaction must be between 0 and 1")

  # if three-way interactions are requested, we need at least 3 base confounders
  if (include_3way && k < 3)
    stop("include_3way requires k >= 3")

  # if two-way interactions are requested, we need at least 2 base confounders
  if (include_2way && k < 2)
    stop("include_2way requires k >= 2")

  # if interaction terms are requested, interaction variance must be > 0
  if ((include_2way || include_3way) && R2_interaction == 0)
    stop("Interaction variance must be > 0 when interaction terms are included")

  # if interaction variance is > 0, at least one interaction toggle must be TRUE
  if (R2_interaction > 0 && !(include_2way || include_3way))
    stop("Interaction variance is positive but no interaction terms are included")

  # number of linear terms (always k main effects)
  k_lin <- k

  # number of interaction terms (two-way and/or three-way, depending on toggles)
  k_int <- 0
  if (include_2way) k_int <- k_int + choose(k, 2)
  if (include_3way) k_int <- k_int + choose(k, 3)

  # split total R^2 into linear and interaction parts
  # of the total variance explained by confounders (R2_total), 1 - R2_interaction is allocated to linear terms
  R2_lin <- (1 - R2_interaction) * R2_total

  # of the total variance explained by confounders (R2_total), R2_interaction is allocated to interaction terms
  R2_int <- R2_interaction * R2_total

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

    # initialize interaction vectors as empty
    dx_int <- dy_int <- numeric(0)

    # sample the delta coefficients for the interaction effects of the confounders on X and Y
    # if we have interaction terms (i.e. k_int > 0) and R2_int > 0, then:
    if (k_int > 0 && R2_int > 0) {

      # sample k_int random numbers from a normal distribution for X
      v <- rnorm(k_int)

      # normalize + scale to hit R2_int
      dx_int <- sqrt(R2_int) * v / sqrt(sum(v^2))

      # sample k_int random numbers from a normal distribution for Y
      v <- rnorm(k_int)

      # normalize + scale to hit R2_int
      dy_int <- sqrt(R2_int) * v / sqrt(sum(v^2))
    }

    # combine linear + interaction parts into full delta vectors
    dx <- c(dx_lin, dx_int)
    dy <- c(dy_lin, dy_int)

    # check if all absolute values are within the specified bounds
    if (all(abs(c(dx, dy)) >= min_abs & abs(c(dx, dy)) <= max_abs)) {

      # if so, return the Delta matrix (2 rows: X and Y)
      D1 <- rbind(dx, dy)

      # set row names
      rownames(D1) <- c("X", "Y")

      # set column names:
      # - first k are main effects: c1..ck
      # - then (if present) interaction names like c1:c2, c1:c3, ..., c1:c2:c3, ...
      if (k_int == 0) {

        # only linear terms
        colnames(D1) <- paste0("c", 1:k)

      } else {

        # build two-way interaction names in the same order as combn(1:k, 2)
        int2_names <- if (include_2way) {
          combn(1:k, 2, FUN = function(ix)
            paste0("c", paste(ix, collapse=":"))
          )
        } else character(0)

        # build three-way interaction names in the same order as combn(1:k, 3)
        int3_names <- if (include_3way) {
          combn(1:k, 3, FUN = function(ix)
            paste0("c", paste(ix, collapse=":"))
          )
        } else character(0)

        # concatenate main effect names + interaction names
        colnames(D1) <- c(paste0("c", 1:k), int2_names, int3_names)
      }

      return(D1)
    }
  }

  # if the loop doesn't return a valid Delta matrix, throw an error
  stop("Failed to sample Delta at t=1 within max_tries")
}