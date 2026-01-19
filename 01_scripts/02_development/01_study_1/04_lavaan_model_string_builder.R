# we want our models to adapt to the number of time points T, and since those models are strings
# we will need to built them using text manipulation
# this is a collection of functions building the following models:
# - CLPM without confounder adjustment
# - CLPM with direct confounder adjustment
# - RI-CLPM with indirect confounder adjustment via random intercepts
# - DPM
# ------------------------------------------------------------------------------------------------------------

# CLPM model string builder, without confounder adjustment at all
build_clpm <- function(T) {

  # here we build the lines:
  # X_t = X_{t-1} + Y_{t-1}
  # Y_t = X_{t-1} + Y_{t-1}
  regress_block <- paste(

    # for each time point from 2 to T
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )

    # add a line break between each time point
    })), collapse="\n"
  )

  # now we need to add the residual covariances
  # producing X_t ~~ Y_t
  resid_cov <- paste(sprintf("x%d ~~ y%d", 1:T, 1:T), collapse="\n")

  # the residual variances for X_t and Y_t
  resid_vars <- paste(

    # yielding lines like X_t ~~ X_t
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),

    # and Y_t ~~ Y_t
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  # we now need to set the means to 1
  means_block <- paste(

    # produces lines: x1 + x2 + ... + xT ~ 1
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",

    # produces lines: y1 + y2 + ... + yT ~ 1
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  # combine all blocks into one model string
  paste(regress_block, resid_cov, resid_vars, means_block, sep="\n")
}

# same as above, but with direct confounder adjustment added
build_clpm_with_Cs <- function(T, k) {

  # creates the line c1 + c2 + ... + ck
  C_names <- paste0("c", 1:k, collapse=" + ")

  # autoregressive and cross-lagged paths, but also confounders added
  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # produces: X_t ~ X_{t-1} + Y_{t-1} + c1 + c2 + ... + ck
        sprintf("x%d ~ x%d + y%d + %s", t, t-1, t-1, C_names),

        # produces: Y_t ~ X_{t-1} + Y_{t-1} + c1 + c2 + ... + ck
        sprintf("y%d ~ x%d + y%d + %s", t, t-1, t-1, C_names)
      )
    })), collapse="\n"
  )

  # from here the function behaves the same as above
  resid_cov <- paste(sprintf("x%d ~~ y%d", 1:T, 1:T), collapse="\n")

  resid_vars <- paste(
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  means_block <- paste(
    paste(paste0("x",1:T), collapse=" + "), "~ 1\n",
    paste(paste0("y",1:T), collapse=" + "), "~ 1\n"
  )

  paste(regress_block, resid_cov, resid_vars, means_block, sep="\n")
}

# same as above, but with indirect confounder adjustment via random intercepts
build_riclpm <- function(T) {

  # here we create the random intercepts
  ri_block <- paste0(

    # produces lines like rix =~ 1*x1 + 1*x2 + ... + 1*xT
    "rix =~ ", paste(sprintf("1*x%d", 1:T), collapse=" + "), "\n",

    # produces lines like riy =~ 1*y1 + 1*y2 + ... + 1*yT
    "riy =~ ", paste(sprintf("1*y%d", 1:T), collapse=" + "), "\n",

    # since this is allways the same, we directly add the variances and covariance of the random intercepts
    "rix ~~ rix\n riy ~~ riy\n rix ~~ riy\n"
  )

  # here we fix the residual variances to zero
  resid_fix <- paste0(

    # produces lines like x1 ~~ 0*x1 + 0*x2 + ... + 0*xT
    paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse="; "), "\n",

    # and y1 ~~ 0*y1 + 0*y2 + ... + 0*yT
    paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the within-person latent variables for X_t and Y_t
  within_lat <- paste0(

    # produces lines like wx1 =~ 1*x1, wx2 =~ 1*x2, ..., wxT =~ 1*xT
    paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse="; "), "\n",

    # and wy1 =~ 1*y1, wy2 =~ 1*y2, ..., wyT =~ 1*yT
    paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the orthogonality constraints: i.e. stable traits are uncorrelated with within-person fluctuations
  orth <- paste0(
    "rix ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "rix ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wx%d", 1:T), collapse=" + "), "\n",
    "riy ~~ ", paste(sprintf("0*wy%d", 1:T), collapse=" + "), "\n"
  )

  # here we create the within-person variances
  within_var <- paste0(

    # creates lines like wx1 ~~ wx1, wx2 ~~ wx2, ..., wxT ~~ wxT
    paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n",

    # and wy1 ~~ wy1, wy2 ~~ wy2, ..., wyT ~~ wyT
    paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the within-person covariances
  within_cov <- paste0(

    # creates lines like wx1 ~~ wy1, wx2 ~~ wy2, ..., wxT ~~ wyT
    paste(sprintf("wy%d ~~ wx%d", 1:T, 1:T), collapse="; "), "\n"
  )

  # here we create the autoregressive and cross-lagged paths
  regress <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}: wx_t ~ wx_{t-1} + wy_{t-1}
        sprintf("wx%d ~ wx%d + wy%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}: wy_t ~ wx_{t-1} + wy_{t-1}
        sprintf("wy%d ~ wx%d + wy%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  # here we create the means
  means <- paste0(

    # produces lines like x1 ~ mx*1, y1 ~ my*1
    paste(paste0("x",1:T), collapse=" + "), " ~ mx*1\n",

    # produces lines like x1 ~ mx*1, y1 ~ my*1
    paste(paste0("y",1:T), collapse=" + "), " ~ my*1\n"
  )

  # finally, we put it all together
  paste(ri_block, resid_fix, within_lat, orth,
        within_var, within_cov, regress, means, sep="\n")
}

# now we build the DPM model string builder
build_dpm <- function(T) {

  # define the accumulating factors FX 
  FX_block <- paste0(

    # produces line FX =~ 1*x1 + 1*x2 + ... + 1*xT
    "FX =~ ", paste(sprintf("1*x%d", 2:T), collapse=" + "), "\n"
  )

  # define the accumulating factors FY
  FY_block <- paste0(

    # produces line FY =~ 1*y1 + 1*y2 + ... + 1*yT
    "FY =~ ", paste(sprintf("1*y%d", 2:T), collapse=" + "), "\n"
  )

  # define the residual covariances between FX and x1, and FY and y1
  fx_cov_block <- "FX ~~ x1 + y1\n"
  fy_cov_block <- "FY ~~ x1 + y1\n"

  # define the autoregressive and cross-lagged paths
  regress_block <- paste(
    unlist(lapply(2:T, function(t){
      c(

        # X_t regressed on X_{t-1} and Y_{t-1}
        sprintf("x%d ~ x%d + y%d", t, t-1, t-1),

        # Y_t regressed on X_{t-1} and Y_{t-1}
        sprintf("y%d ~ x%d + y%d", t, t-1, t-1)
      )
    })), collapse="\n"
  )

  # define the residual covariances between X_t and Y_t
  resid_cov_block <- paste(

    # produces lines like X_t ~~ Y_t
    sprintf("x%d ~~ y%d", 1:T, 1:T),
    collapse="\n"
  )

  # define the latent covariances between FX and FY
  latent_cov_block <- paste(
    "FX ~~ FX",
    "FY ~~ FY",
    "FX ~~ FY",
    sep="\n"
  )

  # define the residual variances
  resid_var_block <- paste(

    # produces lines like X_t ~~ X_t
    paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse="\n"),

    # produces lines like Y_t ~~ Y_t
    paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse="\n"),
    sep="\n"
  )

  # define the means
  means_block <- paste(

    # produces lines like x1 ~ 1, y1 ~ 1
    paste(sprintf("x%d", 1:T), collapse=" + "), "~ 1\n",

    # produces lines like x1 ~ 1, y1 ~ 1
    paste(sprintf("y%d", 1:T), collapse=" + "), "~ 1\n"
  )

  # finally, we put it all together
  paste(
    FX_block,
    FY_block,
    fx_cov_block,
    fy_cov_block,
    regress_block,
    resid_cov_block,
    latent_cov_block,
    resid_var_block,
    means_block,
    sep="\n"
  )
}