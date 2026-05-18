# For BCA (Baseline Covariate Adjustment) in panel data, we need to residualise each observed variable at each time point.
# Concretely, we fit the models: x_t ~ f(U) + e, and y_t ~ f(U) + e, where U are the baseline confounders with linear effects. 
# We then replace the observed x_t with e (the residuals), and similarly for y_t. 
# --------------------------------------------------------------------------------

residualise_panel_linearC <- function(df,
                                      x_prefix = "x",
                                      y_prefix = "y",
                                      c_prefix = "c") {
  
  # convert to data frame
  df <- as.data.frame(df)
  
  # get column names
  x_cols <- grep(paste0("^", x_prefix, "\\d+$"), names(df), value=TRUE)
  y_cols <- grep(paste0("^", y_prefix, "\\d+$"), names(df), value=TRUE)
  c_cols <- grep(paste0("^", c_prefix, "\\d+$"), names(df), value=TRUE)

  # stop if no confounders found
  if (length(c_cols) == 0)
    stop("No confounder columns found.")

  # convert confounders to matrix
  C <- as.matrix(df[c_cols])

  # for each x and y, residualise against confounders
  for (x in x_cols)

    # with the linear model: x_t ~ confounders, and replace the column with the residuals
    df[[x]] <- resid(lm(df[[x]] ~ C))

  # same for y
  for (y in y_cols)
    df[[y]] <- resid(lm(df[[y]] ~ C))

  # return the residualised data frame
  df
}
