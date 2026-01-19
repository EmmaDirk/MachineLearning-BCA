# this script contains functions to generate beta trajectories such that:
# function 1: the effects of baseline confounders remain constant over time
# function 2: the effects of baseline confounders change in a stepwise manner over time such that 
#             the R2 increases from old_R2 to new_R2 at a specified time point
# ---------------------------------------------------------------------

# function (1) to generate B trajectory where B is constant
generate_B_constant <- function(
    B1,                                                           # baseline B-matrix
    T                                                             # number of time points
){

  # create an emtpy list of length T
  B_list <- vector("list", T)

  # set names for each time point
  names(B_list) <- paste0("t", 1:T)

  # fill the list with copies of B1
  for (t in 1:T) {
    B_list[[t]] <- B1
  }
  
  # return the list of B matrices
  return(B_list)
}

# function (2) to generate B trajectory with a hard step
generate_B_stepwise <- function(
    B1,                                                            # baseline B-matrix
    T,                                                             # number of time points
    step_at = floor(T/2) + 1,                                      # when the step starts (default: second half)
    old_R2 = 0.15,                                                 # baseline R2 (before the step)
    new_R2 = 0.40                                                  # higher R2 (after the step)
){

  # scaling factor such that R2 changes from old_R2 to new_R2
  # scaling factor = sqrt(new_R2 / old_R2)
  scale_factor <- sqrt(new_R2 / old_R2)                            # for old_R2=0.15 and new_R2=0.40 -> ~ 1.633

  # build list of B matrices of length T
  B_list <- vector("list", T)

  # set names for each time point
  names(B_list) <- paste0("t", 1:T)

  # fill the list
  for (t in 1:T) {

    # if we are BEFORE the step: keep B exactly equal to B1
    if (t < step_at) {

      # t1, t2, ..., are just baseline
      B_list[[t]] <- B1

    } else {

      # if we are AT the step or AFTER the step: jump to the higher beta matrix
      # this creates the hard step like:
      # {B1, B1, B1, B1*scale_factor, B1*scale_factor}
      B_list[[t]] <- B1 * scale_factor
    }
  }

  return(B_list)
}