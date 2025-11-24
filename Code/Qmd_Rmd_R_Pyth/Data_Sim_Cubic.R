# Here we create a script to create a data simulation for panel data. 
# We simulate X and Y, having autoregressive effects and a lag-1 crossed effect. 
# X and Y are both affected by a number of confounding variables C, which have stable effects
# over time as they represent time-invariant confounding variables, such as IQ. 
# The data simulation model will have the following properties:
#
# 1. Each variable will have an innovation (or unexplained variance) such that the total variance of each variable is 1.
# 2. Each variable will have a mean of 0. 
# 3. The autoregressive effects will be set to 0.25 for both X and Y.
# 4. The crossed lag-1 effect of X on Y will be set to 0.1.
# 5. The confounding variables will have effects of 0.2 on X and 0.15 Y.
# 6. The innovations will have residual correlations of 0.3 at each time-point, but not across time-points.
#
# Let's rephrase what we are doing here. We have got some knowns, because we specify them. They include the
# effects matrix A, the effects matrix B and the variance matrix psi. We also know that the diagonal of S 
# is supposed to be 1. This leaves us with 2 unknowns. We do not know the variance of the innovations, and we
# do not know the covariance between X and Y. Hence all code is actually just a way of solving for these unknowns
# that we need to simulate data. 

# reproducibility
set.seed(427) 

# load required libraries
library(mvtnorm)
library(psych)                                             # for descriptives
library(lavaan)                                            # for fitting CLPM
library(tidyverse)

### CLPM parameters ###
N <- 100000                                                # number of people
T <- 5                                                     # number of time points

ax <- 0.25                                                 # effect of X_(t-1) on X_t
ay <- 0.25                                                 # effect of Y_(t-1) on Y_t
bx <- 0.10                                                 # effect of Y_(t-1) on X_t
by <- 0.10                                                 # effect of X_(t-1) on Y_t
rho <- 0.30                                                # residual correlation of innovations
                                                          
# create the path matrix A from these effects
A <- matrix(c(ax, by, bx, ay), nrow = 2, byrow = TRUE)

### confounder parameters ###
k_lin   <- 3                                               # number of base/linear confounders
k_quad  <- 3                                               # number of "quadratic" confounders
k_cubic <- 3                                               # number of "cubic" confounders

k <- k_lin + k_quad + k_cubic                              # total number of confounders (linear + quadratic + cubic)

# effects of confounders on X and Y
gamma_x <- c(0.30, 0.20, 0.10,                             # linear terms
             0.05, 0.02, 0.00,                             # quadratic terms
             0.02, 0.1, 0.00)                              # cubic terms

gamma_y <- c(0.30, 0.20, 0.10,                             # linear terms
             0.05, 0.02, 0.00,                             # quadratic terms
             0.02, 0.1, 0.00)                              # cubic terms

# create the full path matrix B from these effects (2 x k)
B <- rbind(gamma_x, gamma_y)

# variances of *base* confounders (C1, C2, C3); they are N(0,1) and independent
tau2 <- rep(1, k_lin)                                     

# covariance matrix for the base confounders (3 x 3 identity here)
Psi_base <- diag(tau2)                                     

# generate base confounders
U_lin <- rmvnorm(N, mean = rep(0, k_lin), sigma = Psi_base)

# generate the quadratic and cubic confounders
# we orthogonalize them wrt the linear confounders
# such that no linear dependency remains between them
U_quad <- apply(U_lin, 2, function(col) {

  # raw quadratic transform of each confounder
  z_raw  <- col^2

  # remove the linear part of col (orthogonalize wrt the linear confounder)
  z_perp <- resid(lm(z_raw ~ col))

  # standardize to mean 0 and variance 1
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

# create cubic confounders
U_cubic <- apply(U_lin, 2, function(col) {

  # raw cubic transform of each confounder
  z_raw  <- col^3

  # remove the linear part of col (orthogonalize wrt the linear confounder)
  z_perp <- resid(lm(z_raw ~ col))

  # standardize to mean 0 and variance 1
  as.numeric(scale(z_perp, center = TRUE, scale = TRUE))
})

# combine into full confounder matrix (N x k)
U <- cbind(U_lin, U_quad, U_cubic)

# empirical covariance matrix of all confounder terms (linear, quadratic, cubic)
Psi  <- cov(U)

### varcov matrixes ###
# we want the variance of each variable to be 1
S_target <- diag(2)                                        # since we have two variables, X and Y, we create a 2x2 identity matrix. 
                                                           # This assumes NO covariance between X and Y further than the effects that were specified before.

# the confounders are the only exogenous variables
# and as such, we can calculate the variance that they contribute to each variable (or occasion of the variable) directly
# which is simply the variance of the confounders*effect size.
# Psi contains their varcov matrix, and B contains the effect sizes. 
# we use var(Ax) = A*var(x)*A' -> B*Psi*B'
# Note: here Psi is the full covariance of (C_lin, C_quad, C_cubic).
S_U  <- B %*% Psi %*% t(B)

# now we want to specify a correlation rho between the innovations of X and Y,
# which means the part that is not explained by the confounders or lagged relationships
# still correlates between X and Y. Our job is to find out: given we specify a value rho,
# denoting the correlation between the innovations of X and Y, what must the covariance between X and Y be,
# while accounting for the covariance that is already induced by the confounders and the lagged relationships.
# In one question: Given S_U, A and our desired rho, what must the covariance between X and Y be in order for
# the system to be stationary? 

find_c <- function(A, S_U, rho) {

  # function of c: difference between actual innovation correlation and target rho
  f <- function(c) {

    # target stationary covariance of (X_t, Y_t)
    # the variances are 1 since we specified them to be
    # and the covariance is c, which we are searching for
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    
    # dynamic part variance
    # first we get rid of the confounder-induced variance
    S_dyn_c <- S_target_c - S_U
    
    # now we can compute the innovation covariance matrix
    # based on the stationarity constraints
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    
    # compute correlation of innovations
    # the variance of innovations of X is v1
    v1 <- Sigma_e_c[1, 1]
    # the variance of innovations of Y is v2
    v2 <- Sigma_e_c[2, 2]
    # their covariance is cov12
    cov12 <- Sigma_e_c[1, 2]
    # their correlation is then:
    corr_e <- cov12 / sqrt(v1 * v2)
    
    # this is then the difference between the desired correlation and the actual correlation
    corr_e - rho
  }

  # now this function f(c) should be zero when we have found the correct c
  # as such we can look for the value of c that makes f(c) = 0 (i.e., corr_e = rho).
  uniroot(f, interval = c(-0.99, 0.99))$root
}

# use the function to find the required covariance between X and Y
c_stat <- find_c(A, S_U, rho)

# overwrite S_target to include the stationary covariance between X and Y
S_target <- matrix(c(1, c_stat,
                     c_stat, 1),
                   nrow = 2, byrow = TRUE)

# variance by the dynamic process:
# now we know that any variance that does not come directly from the confounder, must come from the dynamic process (feed forward + innovations)
# since we also know that we want the variance of our variables to be 1, we can calculate the variance that we can still 'spend'
# on the dynamic part by simple subtraction:
S_dyn <- S_target - S_U

# now we can check if our matrix is positive semidefinite (i.e. that all eigenvalues are positive)
# i.e. we do not want to set our requirement so that we must generate data with negative variances. 
# we check that here.
# first we make it robust to small numerical errors created by floating point storage (S_dyn + t(S_dyn))/2
# and then extract the eigenvalues
eS <- eigen((S_dyn + t(S_dyn))/2)$values
# check if any are very close to, or smaller than 0. 
eS[eS < 1e-10]

# now we need to compute the innovation covariance matrix.
# we compute how much innovation variance we must add, such that the dynamic process is stationary. 
# which means that the variance of X_t and Y_t do not change over time.
# read it like: how much innovation variance must we add such that the variance of each variable equals S_dyn. 
# and remember that if the variance of the dynamic process equals S_dyn, then all our variables have variances equal to S_target.
# 1. X_t = A X_(t-1) + e_t
# 2. for all variables: D_t = A D_{t-1} + e_t
# 3. S_dyn = var(D_t) = A S_dyn A' + S_e (Since var(X_t) = A^2 * var(X_{t-1}) + var(e_t))
# 4. solve for S_e: S_e = S_dyn - A S_dyn A'
# note that this only works because the system is stationary! we can exclude the variance D_t receives directly from the confounders
# since it is the same at every step, and the variance of the dynamic process is the same at every step. If this would not be the case
# we would need to use formulas to compute how much innovation should be added to each measurement occasion separately.
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A

# and again we fix potential rounding error
Sigma_e <- (Sigma_e + t(Sigma_e)) / 2

# For t = 1, there is no carryover from a previous wave. The observed var-cov should be:
# Var([X1,Y1]) = S_U + Sigma_e1. To keep the same S_target at t=1, set:
Sigma_e1 <- S_dyn                             # because S_target = S_U + S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1)) / 2

# now we can simulate the dynamic process data for X and Y
# we add here the entire variance to each variable, so that is feed forward + innovations
# these are only X1 and Y1, so they only have the direct confounder effects

#  for t = 1: no lag yet, so Ddyn = e1 with Sigma_e1 ----
Ddyn <- rmvnorm(N,                            
  mean  = c(0, 0),                            # mean of the variables X and Y = 0, since we want standardised effect sizes
  sigma = Sigma_e1)                           # dynamic part covariance for the first wave

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  paste0("c", 1:k))

# save the confounder
df[, (2*T + 1):(2*T + k)] <- U

# we now add to the first wave X and Y the direct confounder effects
# where we have our confounders (U) and the effect sizes (B)
obs1 <- Ddyn + U %*% t(B)

# we can store the first wave in our dataframe now
df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# a next wave can be simulated as:
# Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma=Sigma_e)
# where we get the previous wave observations, times their effect sizes + the innovations such that var=1. 
# and afterwards we can again add the direct confounder effects
# obs <- Ddyn + U %*% t(B)
# write a loop to do this up to T occasions:
for(i in 2:T){ 

  # simulate the dynamic part
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 

  # add the direct confounder effects
  obs <- Ddyn + U %*% t(B) 

  # move the data to the dataframe
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2]
}

# check that the variance of each variable are 1
round(apply(df, 2, var), 3)

# check the varcov matrix
sigma <- cov(df)
round(sigma, 3)

# check descriptives
describe(df)

# create a sample dataframe for plotting residuals 
df_sample <- df[sample(1:nrow(df), 10000), ]

# fit a linear regression between X1 and C1
lm_fit <- lm(x1 ~ c1, data = as.data.frame(df_sample))

# save the residuals
residuals <- lm_fit$residuals

# plot the residuals versus C1
ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  theme_minimal()

# we will run this script various times using different effect sizes to produce
# multiple datasets with different magnitudes of non-linearity. 
write.csv(df, file = "Thesis/Code/Data/Data_Sim_Cubic_Pop.csv", row.names = FALSE)
