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
# lets start with defining some effect sizes as we do in our scripts

### confounder parameters ###
k_lin   <- 3                                               # number of base/linear confounders
k_quad  <- 3                                               # number of "quadratic" confounders
k_cubic <- 3                                               # number of "cubic" confounders

k <- k_lin + k_quad + k_cubic                              # total number of confounders

# effects of confounders on X and Y
gamma_x <- c(0.30, 0.20, 0.10,                             # linear terms
             0.10, 0.05, 0.00,                             # quadratic terms
             0.05, 0.02, 0.00)                             # cubic terms

gamma_y <- c(0.30, 0.20, 0.10,                             # linear terms
             0.10, 0.05, 0.00,                             # quadratic terms
             0.05, 0.02, 0.00)                             # cubic terms

# now save the linear effects seperately for later
gamma_x_lin <- gamma_x[1:k_lin]                            # linear effects on X
gamma_y_lin <- gamma_y[1:k_lin]                            # linear effects on Y

# and the non-linear effects
gamma_x_nl <- gamma_x[(k_lin + 1):k]                       # non-linear effects on X
gamma_y_nl <- gamma_y[(k_lin + 1):k]                       # non-linear effects on Y

# now write the function that computes the non-linearity measure
nonlinearity_eta <- function(gamma_lin, gamma_nl) {
  var_L <- sum(gamma_lin^2)
  var_N <- sum(gamma_nl^2)
  var_N / (var_L + var_N)
}

# compute non-linearity measures for X and Y
nonlin_X <- nonlinearity_eta(gamma_x_lin, gamma_x_nl)
nonlin_Y <- nonlinearity_eta(gamma_y_lin, gamma_y_nl)

# take their average
mean(c(nonlin_X, nonlin_Y))

# now interpret this as follows
# eta close to 0 means mostly linear
# eta close to 1 means mostly non-linear
# eta over 0.5 means our data follow a really really non-linear pattern

# we can upgrade our function a bit to do all in one step (provided k and k_lin are known)
# so we create the vectors, then compute the non-linearity measure, do this for both X and Y,
# and finally take the average: 

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

# a little more non-linear
gamma_x2 <- c(0.25, 0.15, 0.10,
               0.10, 0.07, 0.05,
               0.05, 0.03, 0.02)

gamma_y2 <- c(0.25, 0.15, 0.10,
               0.10, 0.07, 0.05,
               0.05, 0.03, 0.02)

nonlinearity_index(gamma_x2, gamma_y2, k_lin)  # expect high value

# really quite non-linear
gamma_x3 <- c(0.20, 0.15, 0.10,
               0.15, 0.10, 0.05,
               0.10, 0.07, 0.05)

gamma_y3 <- c(0.20, 0.15, 0.10,
               0.15, 0.10, 0.05,
               0.10, 0.07, 0.05)

nonlinearity_index(gamma_x3, gamma_y3, k_lin)  # expect high value

# extremely non-linear
gamma_x4 <- c(0.10, 0.05, 0.02,
               0.20, 0.15, 0.10,
               0.15, 0.10, 0.05)

gamma_y4 <- c(0.10, 0.05, 0.02,
               0.20, 0.15, 0.10,
               0.15, 0.10, 0.05)

nonlinearity_index(gamma_x4, gamma_y4, k_lin)  # expect very high value
