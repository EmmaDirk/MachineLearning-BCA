# this script contains functions to generate delta trajectories from our baseline sampled Delta-Matrix D1
# these delta trajectories represent how the effects of baseline confounders change over time
# we have two scenarios implemented here:
# - the effects remain constant over time (do not change)
# - the effects change in a stepwise manner at a specified time point
#
# function 1: the effects of baseline confounders remain constant over time
# function 2: the effects of baseline confounders change in a stepwise manner over time such that 
#             the variance explained (R2) by the confounders increases from old_R2 to new_R2 at a 
#             specified time point
# ---------------------------------------------------------------------

# function (1) to generate Delta trajectory where the coefficients are constant
generate_D_constant <- function(
    D1,                                                           # baseline D-matrix
    T                                                             # number of time points
){

  # create an emtpy list of length T
  D_list <- vector("list", T)

  # set names for each time point
  names(D_list) <- paste0("t", 1:T)

  # fill the list with copies of D1
  for (t in 1:T) {
    D_list[[t]] <- D1
  }
  
  # return the list of D matrices
  return(D_list)
}

# function (2) to generate D trajectory with a hard step
generate_D_stepwise <- function(
    D1,                                                            # baseline D-matrix
    T,                                                             # number of time points
    step_at = floor(T/2) + 1,                                      # when the step starts (default: second half)
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40                                                  # higher (or lower) R2 (after the step)
){

  # scaling factor such that R2 changes from old_R2 to new_R2
  # scaling factor = sqrt(new_R2 / old_R2)
  scale_factor <- sqrt(new_R2 / old_R2)                            # for old_R2=0.15 and new_R2=0.40 -> ~ 1.633

  # build list of D matrices of length T
  D_list <- vector("list", T)

  # set names for each time point
  names(D_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep D exactly equal to D1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      D_list[[t]] <- D1

    } else {

      # if we are AT the step or AFTER the step: jump to the higher delta matrix
      # this creates the hard step like:
      # {D1, D1, D1, D1*scale_factor, D1*scale_factor}
      D_list[[t]] <- D1 * scale_factor
    }
  }

  # return the list of D matrices
  return(D_list)
}
