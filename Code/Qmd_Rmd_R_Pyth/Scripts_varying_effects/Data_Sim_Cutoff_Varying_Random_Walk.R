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
k_lin <- 3                                                 # number of "linear" confounders (baseline C1..C3)

# baseline linear effects of confounders on X_t and Y_t
gamma_x_lin <- c(0.30, 0.20, 0.10)                         # effect of confounders on X_t
gamma_y_lin <- c(0.30, 0.20, 0.10)                         # effect of confounders on Y_t

# cutoff for plateau effect
thr <- 1                                                   # cutoff threshold for confounders (in SD units)

# effect of the "excess above the threshold" on X and Y
gamma_x_cut <- -gamma_x_lin                                # plateau for X above threshold
gamma_y_cut <- -gamma_y_lin                                # plateau for Y above threshold

# total number of confounder variables (C and cutoff versions)
k <- k_lin * 2                                             

# simulate exogenous linear confounders
Psi_base <- diag(rep(1, k_lin))                            
U_lin <- rmvnorm(N, 
                 mean  = rep(0, k_lin),                    
                 sigma = Psi_base)                         

colnames(U_lin) <- paste0("c", 1:k_lin)                    # C1, C2, C3

# cutoff confounders: H_j = max(C_j - thr, 0)
U_cut <- pmax(U_lin - thr, 0)
colnames(U_cut) <- paste0("c", 1:k_lin, "_cut")

# combine linear and cutoff confounders into one matrix U
U <- cbind(U_lin, U_cut)                                   # N x k, where k = 2 * k_lin

# empirical varcov matrix
Psi <- cov(U)

# effect vectors including cutoff effects
# order: (C1, C2, C3, H1, H2, H3)
gamma_x <- c(gamma_x_lin, gamma_x_cut)
gamma_y <- c(gamma_y_lin, gamma_y_cut)

### RANDOM WALK FOR BETA_t WITH SHRINKAGE TOWARD B0 ###

# baseline 2 x k effect matrix from all confounders to (X, Y)
B0 <- rbind(gamma_x, gamma_y)                              # 2 x k

# RW mean vector for each coefficient (trend); here zero
RW_mean <- matrix(0, nrow = 2, ncol = k)

# random walk sd and shrinkage toward B0
rw_sd  <- 0.03                                             # step variability
shrink <- 0.10                                             # pull toward B0

Beta_list <- vector("list", T)
Beta_list[[1]] <- B0

# create random walk with shrinkage toward B0
for(t in 2:T){

  # random walk step with mean vector
  step <- matrix(
    rnorm(length(B0), mean = as.vector(RW_mean), sd = rw_sd),
    nrow = 2
  )

  B_prev <- Beta_list[[t-1]]

  # RW step + shrinkage toward B0
  B_new <- B_prev + step + shrink * (B0 - B_prev)

  Beta_list[[t]] <- B_new
}

### varcov matrices ###

# want variance of each variable to be 1
S_target <- diag(2)

# variance contributed by confounders over time:
# S_U_t = B_t * Psi * B_t'
S_U_list <- lapply(Beta_list, function(Bt) Bt %*% Psi %*% t(Bt))
S_U  <- Reduce("+", S_U_list) / length(S_U_list)           # mean across time

# find covariance c between X and Y needed for target residual correlation rho
find_c <- function(A, S_U, rho) {

  f <- function(c) {

    # target stationary covariance of (X_t, Y_t)
    S_target_c <- matrix(c(1, c,
                           c, 1),
                         nrow = 2, byrow = TRUE)
    
    # dynamic part variance
    S_dyn_c <- S_target_c - S_U
    
    # innovation covariance matrix under stationarity
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    
    # correlation of innovations
    v1 <- Sigma_e_c[1, 1]
    v2 <- Sigma_e_c[2, 2]
    cov12 <- Sigma_e_c[1, 2]
    corr_e <- cov12 / sqrt(v1 * v2)
    
    corr_e - rho
  }

  uniroot(f, interval = c(-0.99, 0.99))$root
}

# required covariance between X and Y
c_stat <- find_c(A, S_U, rho)

# stationary covariance of (X, Y)
S_target <- matrix(c(1, c_stat,
                     c_stat, 1),
                   nrow = 2, byrow = TRUE)

# variance by the dynamic process:
S_dyn <- S_target - S_U
S_dyn <- (S_dyn + t(S_dyn))/2                             # symmetrise

# check positive semidefinite
eS <- eigen(S_dyn)$values
eS[eS < 1e-10]

# innovation covariance matrix
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2

# For t = 1, Var([X1,Y1]) = S_U + Sigma_e1 = S_target
Sigma_e1 <- S_dyn
Sigma_e1 <- (Sigma_e1 + t(Sigma_e1))/2

# simulate dynamic process

# t = 1: no lag yet, so Ddyn = e1 with Sigma_e1
Ddyn <- rmvnorm(N,                            
  mean  = c(0,0),                              
  sigma = Sigma_e1)                            

# initialise df with N rows and 2*T + k columns. 2*T for X and Y, k for confounders
df <- matrix(NA, nrow = N, ncol = 2*T + k)
colnames(df) <- c(paste0("x", 1:T),
                  paste0("y", 1:T),
                  colnames(U))

# save the confounders (C1..C3 and H1..H3)
df[, (2*T + 1):(2*T + k)] <- U

# first wave: add direct confounder effects using Beta_list[[1]]
obs1 <- Ddyn + U %*% t(Beta_list[[1]])

df[, 1]    <- obs1[, 1] 
df[, 1+T ] <- obs1[, 2]

# remaining waves: time-varying Beta_list[[i]]
for(i in 2:T){ 

  # simulate the dynamic part
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma = Sigma_e) 

  # add the direct confounder effects
  obs <- Ddyn + U %*% t(Beta_list[[i]]) 

  # move the data to the dataframe
  df[, i]    <- obs[, 1] 
  df[, i+T ] <- obs[, 2] 
}

df <- as.data.frame(df)

# checks
round(apply(df, 2, var), 3)
sigma <- cov(df)
round(sigma, 3)
describe(df)

# residual nonlinearity check
df_sample <- df[sample(1:nrow(df), 10000), ]
lm_fit    <- lm(x1 ~ c1, data = df_sample)
residuals <- lm_fit$residuals

ggplot(df_sample, aes(x = c1, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# cache data
write.csv(df, file = "Thesis/Code/Data/Data_Sim_plateau_RW_Pop.csv", row.names = FALSE)
