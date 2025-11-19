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
A <- matrix(c(ax, by, 
              bx, ay), 
            nrow = 2, byrow = TRUE)

### confounder parameters ###
k_lin <- 5                                                 # number of "basic" linear confounders (C1..C5)

# effects of linear confounders on X_t and Y_t
gamma_x_lin <- c(0.30, 0.25, 0.20, 0.15, 0.10)             # effect of confounders on X_t
gamma_y_lin <- c(0.30, 0.25, 0.20, 0.15, 0.10)             # effect of confounders on Y_t

# next we want to add a form of non-linearity. Namely three-way interactions between the confounders.
# We will create new confounders C_ijk = Ci * Cj * Ck for all 3-way combinations of the 5 confounders.
# This is still a "linear" effect in the model, but in the original confounders it is a non-linear term.

# total number of three-way interaction terms for 5 confounders is "5 choose 3" = 10
k_int <- choose(k_lin, 3)                                  # number of three-way interaction terms
k <- k_lin + k_int                                         # total number of confounder variables

# we specify the effect of each three-way interaction on X and Y
gamma_x_int <- rep(0.10, k_int)                            # effect of each 3-way interaction on X_t
gamma_y_int <- rep(0.10, k_int)                            # effect of each 3-way interaction on Y_t

# combine all effects into gamma vectors
# order will be: (C1, C2, C3, C4, C5, then all C_ijk interactions)
gamma_x <- c(gamma_x_lin, gamma_x_int)
gamma_y <- c(gamma_y_lin, gamma_y_int)

# path matrix B from confounders to (X, Y), 2 x k
B <- rbind(gamma_x, gamma_y)

# we can first generate our exogenous linear confounders, since they are exogenous
# we use the function rmvnorm to generate multivariate normal data
Psi_base <- diag(rep(1, k_lin))                            # variance-covariance matrix of base confounders (I in this case)
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    # mean of confounders = 0, since we want standardised effect sizes
                 sigma = Psi_base)                         # base confounders C1..C5 ~ N(0,1), independent

# now we create the three-way interaction terms:
# for all triples (i, j, k) of the 5 confounders:
comb_3 <- combn(k_lin, 3)                                  # 3 x k_int matrix of indices
U_int_raw <- matrix(NA, nrow = N, ncol = k_int)
for (idx in 1:k_int) {
  i <- comb_3[1, idx]
  j <- comb_3[2, idx]
  k3 <- comb_3[3, idx]
  # C_ijk = Ci * Cj * Ck
  U_int_raw[, idx] <- U_lin[, i] * U_lin[, j] * U_lin[, k3]
}

# we can standardise the interaction terms so that they have mean 0 and variance 1
# (this makes the interpretation of gamma_x_int, gamma_y_int as "standardised" effect sizes easier)
U_int <- scale(U_int_raw)

# now we can combine the linear confounders and the interaction terms into one matrix U
# columns: first the k_lin linear confounders, then the k_int three-way interactions
U <- cbind(U_lin, U_int)                                   # N x k

# and empirically save their varcov matrix
Psi <- cov(U)                                              # full covariance of (C1..C5, all 3-way interactions)

### varcov matrixes ###
# we want the variance of each variable to be 1
S_target <- diag(2)                                        # since we have two variables, X and Y, we create a 2x2 identity matrix. This assumes NO covariance between X and Y
                                                           # further than the effects that were specified before

# the confounders are the only exogenous variables
# and as such, we can calculate the variance that they contribute to each variable (or occasion of the variable) directly
# which is simply the variance of the confounders*effect size 
# Psi contains their varcov matrix, and B contains the effect sizes. 
# we use var(Ax) = A*var(x)*A' -> B*Psi*B'.
# Note that Psi != I here, since the interactions are deterministic non-linear functions of C1..C5,
# which creates covariance between the linear terms and the interaction terms.
S_U  <- B %*% Psi %*% t(B)

# now we want to specify a correlation rho between the innovations of X and Y.
# which means the part that is not explained by the confounders or lagged relationships
# still correlate with each other. Our job is to find out: given we specify a value rho,
# denoting the correlation between the innovations of X and Y, what must the covariance between X and Y be,
# while accounting for the covariance that is already induced by the confounders and the lagged relationships.
# In one question: Given S_U, A and our desired rho, what must the covariance between X and Y be in order for
# the system to be stationary? 

find_c <- function(A, S_U, rho) {

  # function of c: difference between actual innovation correlation and target rho
  # we need to calculate how much covariance is induced by c and the lagged relationships
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
  # as such we can look for the value of c that makes f(c) = 0.
  # Note that if the matrix A would be symmetric, then v1=v2, and we have
  # a closed form solution for v. With two unknowns, we lack an equation and
  # need to solve numerically. 
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
eS[eS<1e-10]

# now we need to compute the innovation covariance matrix.
# we compute how much innovation variance must be added, such that the dynamic process is stationary. 
# which means that the variance of X_t and Y_t do not change over time.
# read it like: how much innovation variance must we add such that the variance of each variable equals S_dyn. 
# and remember that if the variance of the dynamic process equals S_dyn, then all our variables have variances equal to S_target.
# 1. X_t = A X_(t-1) + e_t
# 2. for all variables: D_t = A D_t-1 + e_t
# 3. S_dyn = var(D_t) = A S_dyn A' + S_e (Since var(X_t) = A^2 * var(X_t-1) + var(e_t))
# 4. solve for S_e: S_e = S_dyn - A S_dyn A'
# note that this only works because the system is stationary! we can exclude the variance D_t receives directly from the confounders
# since it is the same at every step, and the variance of the dynamic process is the same at every step. If this would not be the case
# we would need to use formulas to compute how much innovation should be added to each measurement occasion separately.
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A

# and again we fix potential rounding error
Sigma_e <- (Sigma_e + t(Sigma_e))/2

# For t = 1, there is no carryover from a previous wave. The observed var-cov should be:
# Var([X1,Y1]) = S_U + Sigma_e1. To keep the same S_target at t=1, set:
Sigma_e1 <- S_dyn                             # because S_target = S_U + S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

# now we can simulate the dynamic process data for X and Y
# we add here the entire variance to each variable, so that is feed forward + innovations
# these are only X1 and Y1, so they only have the direct confounder effects

#  for t = 1: no lag yet, so Ddyn = e1 with Sigma_e1 ----
Ddyn <- rmvnorm(N,                            
  mean  = c(0,0),                              # mean of the variables X and Y = 0, since we want standardised effect sizes
  sigma = Sigma_e1)                            # dynamic part covariance for the first wave

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  paste0("c", 1:k))

# save the confounders
# c1..c5 are the linear Cs, c6..c15 are the three-way interactions
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

# example: fit a linear regression between X1 and C1
lm_fit <- lm(x1 ~ c1, data = as.data.frame(df_sample))

# save the residuals
residuals <- lm_fit$residuals

# plot the residuals versus C1
ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# we will run this script various times using different effect sizes to produce
# multiple datasets with different magnitudes of non-linearity. 
# write.csv(df, file = "Code/Data/Data_Sim_threeway5_Pop.csv", row.names = FALSE)
