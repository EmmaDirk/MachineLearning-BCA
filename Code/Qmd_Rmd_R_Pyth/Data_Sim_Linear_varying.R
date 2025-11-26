# Here we create a script to create a data simulation for panel data. 
# We simulate X and Y, having autoregressive effects and a lag-1 crossed effect. 
# X and Y are both affected by a number of confounding variables C, which have stable effects
# over time as they represent time-invariant confounding variables, such as IQ. 
# The data simulation model will have the following properties:

# 1. Each variable will have an innovation (or unexplained variance) such that the total variance of each variable is 1.
# 2. Each variable will have a mean of 0. 
# 3. The autoregressive effects will be set to 0.25 for both X and Y.
# 4. The crossed lag-1 effect of X on Y will be set to 0.1.
# 5. The confounding variables will have effects of 0.2 on  X and 0.15 Y.
# 6. The innovations will have residual correlations of 0.3 at each time-point, but not across time-points.

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
A <- matrix(c(ax, by, bx, ay), nrow=2, byrow=TRUE)

### confounder parameters ###
k <- 3                                                     # number of confounders
gamma_x <- c(0.30, 0.20, 0.10)                             # effect of confounders on X_t, where the first entry is the effect of confounder 1 on X_t, etc. 
gamma_y <- c(0.30, 0.20, 0.10)                             # effect of confounders on Y_t, where the first entry is the effect of confounder 1 on Y_t, etc.
tau2 <- c(1, 1, 1)                                         # variances of confounders

B <- rbind(gamma_x, gamma_y)                               # create the path matrix B from these effects (baseline, wave 1)
Psi  <- diag(tau2)                                         # variance-covariance matrix of confounders, which is just a k*k identity matrix since var=1 and cov=0

# now we want to vary these confounder parameters over the time points, to create a time-invariant confounder with varying effects
# we want to do this in a parameterised way, so we can easily change the effects later on. As such we will pick the change of a time-point
# from the previous time-point from a normal distribution with a variance and mean. The variance determines how unstable the effect of the 
# confounder is and the mean determines the direction of change. The mean vector can be fully specified, while the variance is the same
# for all time points.

# IMPORTANT: these vectors must have length T. Entry t is the mean change that is applied when going from wave (t-1) to wave t.
# Example below: a "COVID-like" shock at wave 3 (large negative change), followed by recovery.
mean_change_x_vec <- c(0.00, 0.00, -0.30, 0.05, 0.05)      # change of the mean of the normal distribution from which we draw the change in effect sizes for X
mean_change_y_vec <- c(0.00, 0.00, -0.20, 0.03, 0.03)      # change of the mean of the normal distribution from which we draw the change in effect sizes for Y

sd_change     <- 0.005                                     # standard deviation of the normal distribution from which we draw the change in effect sizes
                                                           # larger values -> more unstable/confounded effects over time

# list to store B per wave
B_list <- vector("list", T)

# first wave is just the original effects
B_list[[1]] <- rbind(gamma_x, gamma_y)

# create time-varying effects
# we build the effects over time as a random walk:
# B_t = B_(t-1) + Delta_t, where Delta_t ~ Normal(mean_change_vec[t], sd_change^2)
for(t in 2:T){

  # extract the mean changes for this time point
  mean_change_x <- mean_change_x_vec[t]
  mean_change_y <- mean_change_y_vec[t]

  # draw random changes for each coefficient of X and Y at time t
  change_x <- rnorm(k, mean = mean_change_x, sd = sd_change)
  change_y <- rnorm(k, mean = mean_change_y, sd = sd_change)

  # apply changes to previous wave's coefficients
  new_gamma_x <- B_list[[t-1]][1, ] + change_x
  new_gamma_y <- B_list[[t-1]][2, ] + change_y

  # save as 2xk matrix for wave t
  B_list[[t]] <- rbind(new_gamma_x, new_gamma_y)
}

### helper function to find covariance c, given S_U and target innovation correlation rho ###
# now we want to specify a correlation rho between the innovations of X and Y.
# which means the part that is not explained by the confounders or lagged relationships
# still correlate with each other. Our job is to find out: given we specify a value rho,
# denoting the correlation between the innovations of X and Y, what must the covariance between X and Y be,
# while accounting for the covariance that is already induced by the confounders and the lagged relationships.
# In one question: Given S_U, A and our desired rho, what must the covariance between X and Y be in order for
# the system to be stationary (for that given S_U)? 

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
    # the variation of innovations of X is v1
    v1 <- Sigma_e_c[1, 1]
    # the variation of innovations of Y is v2
    v2 <- Sigma_e_c[2, 2]
    # their covariance is cov12
    cov12 <- Sigma_e_c[1, 2]
    #their correlation is then:
    corr_e <- cov12 / sqrt(v1 * v2)
    
    # this is then the difference between the desired correlation and the actual correlation
    corr_e - rho
  }

  # now this function f(c) should be zero when we have found the correct c
  # as such we can look for the value of corr_e that makes f(c) = 0.
  # Note that if the matrix A would be symmetric, then v1=v2, and we have
  # a closed form solution for v. With two uknowns, we lack an equation and
  # need to solve numerically. 
  uniroot(f, interval = c(-0.99, 0.99))$root
}

### time-point specific variance-covariance matrices ###

# here we enforce that at EACH time point:
#  1) Var(X_t) = Var(Y_t) = 1
#  2) Corr(innovation_X_t, innovation_Y_t) = rho
# given:
#  - time-varying confounder effects B_list[[t]]
#  - autoregressive and cross-lagged structure in A

# we will store, for each time point:
#  - S_target_list[[t]] : target var-cov matrix of (X_t, Y_t)
#  - S_U_list[[t]]      : variance induced by confounders at time t
#  - S_dyn_list[[t]]    : dynamic part variance at time t
#  - Sigma_e_list[[t]]  : innovation covariance at time t
#  - Sigma_e1           : innovation covariance at t=1 (no lag, only innovation)

S_target_list <- vector("list", T)
S_U_list      <- vector("list", T)
S_dyn_list    <- vector("list", T)
Sigma_e_list  <- vector("list", T)

# loop over time points to compute all components
for(t in 1:T){

  # confounder-induced variance at time t
  B_t   <- B_list[[t]]
  S_U_t <- B_t %*% Psi %*% t(B_t)
  S_U_list[[t]] <- S_U_t
  
  # find the stationary covariance c_t between X and Y at time t
  # such that the residual (innovation) correlation equals rho, 
  # given S_U_t and A.
  c_t <- find_c(A, S_U_t, rho)
  
  # target stationary covariance matrix at time t
  S_target_t <- matrix(c(1, c_t,
                         c_t, 1),
                       nrow = 2, byrow = TRUE)
  S_target_list[[t]] <- S_target_t
  
  # dynamic part variance at time t
  S_dyn_t <- S_target_t - S_U_t
  # symmetrise for numerical stability
  S_dyn_t <- (S_dyn_t + t(S_dyn_t))/2
  S_dyn_list[[t]] <- S_dyn_t
  
  # innovation covariance at time t
  Sigma_e_t <- S_dyn_t - t(A) %*% S_dyn_t %*% A
  # fix potential rounding error
  Sigma_e_t <- (Sigma_e_t + t(Sigma_e_t))/2
  Sigma_e_list[[t]] <- Sigma_e_t
}

# For t = 1, there is no carryover from a previous wave. 
# The observed var-cov should be:
# Var([X1,Y1]) = S_U1 + Sigma_e1. To keep S_target_1 at t=1, set:
Sigma_e1 <- S_dyn_list[[1]]                 # because S_target_1 = S_U1 + S_dyn_1
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2      # symmetrise

# For t >= 2, Sigma_e_list[[t]] is already computed above and will be used in the dynamic recursion.

# we can now check that our S_dyn_t matrices are positive semidefinite
eS_all <- lapply(S_dyn_list, function(M){
  eigen((M + t(M))/2)$values
})
# inspect if any eigenvalues are very small or negative
unlist(eS_all)[unlist(eS_all) < 1e-10]

# we can first generate our confounders, since they are exogenous
# we use the function rmvnorm to generate multivariate normal data
U <- rmvnorm(N, 
  mean=rep(0, k),                             # mean of confounders = 0, since we want standardised effect sizes
  sigma=Psi)                                  # variance-covariance matrix of confounders (I in this case)

# now we can simulate the dynamic process data for X and Y
# we add here the entire variance to each variable, so that is feed forward + innovations
# these are only X1 and Y1, so they only have the direct confounder effects

#  for t = 1: no lag yet, so Ddyn = e1 with Sigma_e1 ----
Ddyn <- rmvnorm(N,                            
  mean=c(0,0),                                # mean of the variables X and Y = 0, since we want standardised effect sizes
  sigma=Sigma_e1)                             # dynamic part covariance for the first wave

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  paste0("c", 1:k))

# save the confounder
df[, (2*T + 1):(2*T + k)] <- U

# we now add to the first wave X and Y the direct confounder effects
# where we have our confounders (U) and the effect sizes (B_list[[1]])
obs1 <- Ddyn + U %*% t(B_list[[1]])

# we can store the first wave in our dataframe now
df[, 1] <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# a next wave can be simulated as:
# Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma=Sigma_e_list[[i]])
# where we get the previous wave observations, times their effect sizes + the innovations such that var=1. 
# and afterwards we can again add the direct confounder effects
# obs <- Ddyn + U %*% t(B_list[[i]])
# write a loop to do this up to T occasions:
for(i in 2:T){ 

  # simulate the dynamic part for wave i
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e_list[[i]]) 

  # add the direct confounder effects for time point i (time-varying B)
  obs <- Ddyn + U %*% t(B_list[[i]]) 

  # move the data to the dataframe
  df[, i]   <- obs[, 1] 
  df[, i+T] <- obs[, 2] 
}

# check that the variance of each variable are 1 (approximately, up to numerical error)
round(apply(df, 2, var), 3)

# check the varcov matrix
sigma <- cov(df)
round(sigma, 3)

# check descriptives
describe(df)

# create a sample dataframe for plotting residuals 
df_sample <- df[sample(1:nrow(df), 10000), ]

# fit a linear regression between X at t=1 and Y at t=2
lm_fit <- lm(x1 ~ c1, data = as.data.frame(df_sample))

# save the residuals
residuals <- lm_fit$residuals

# plot the residuals versus X1
ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# cache the data to the data folder
# write.csv(df, file = "Thesis/Code/Data/Data_Sim_linear_Pop.csv", row.names = FALSE)
