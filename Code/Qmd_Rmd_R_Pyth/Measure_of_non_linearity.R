# it helps to have some measure of non-linearity of a relationship
# in our case where we simulate non-linearities always by adding extra variables
# we can measure the ammount of non-linearity directly by looking at the effect sizes
# of the non-linear variables in relation to their linear counterparts

# we say:
# X_t = L(C) + N(C) + E

# we can define a crude non-linearity measure by taking:
# eta_X = var(L(C)) / var(L(C) + N(C))
# then we do the same for Y and average them to a scalar value
# we can make this trivially simple since the variance of each component is 1

nonlinearity_index <- function(gamma_x, gamma_y, k_lin) {
  
  # helper
  eta <- function(g_lin, g_nl) {
    var_L <- sum(g_lin^2)
    var_N <- sum(g_nl^2)
    var_N / (var_L + var_N)
  }
  
  # split into linear and nonlinear parts
  k_total <- length(gamma_x)
  gamma_x_lin <- gamma_x[1:k_lin]
  gamma_y_lin <- gamma_y[1:k_lin]
  
  gamma_x_nl  <- gamma_x[(k_lin + 1):k_total]
  gamma_y_nl  <- gamma_y[(k_lin + 1):k_total]
  
  # compute nonlinearity for X and Y
  eta_X <- eta(gamma_x_lin, gamma_x_nl)
  eta_Y <- eta(gamma_y_lin, gamma_y_nl)
  
  # return their average
  mean(c(eta_X, eta_Y))
}

# test it
nonlinearity_index(gamma_x, gamma_y, k_lin)

# lets now produce a few different gamma matrixes and test the function
# mostly linear
gamma_x1 <- c(0.30, 0.20, 0.10,
               0.05, 0.02, 0.00,
               0.02, 0.01, 0.00)

gamma_y1 <- c(0.30, 0.20, 0.10,
               0.05, 0.02, 0.00,
               0.02, 0.01, 0.00)

nonlinearity_index(gamma_x1, gamma_y1, k_lin)  # expect low value

