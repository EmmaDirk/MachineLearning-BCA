# this script extends the previous data simulation script by including effect from confounders that vary over time
# do not that we still want a stationary system and the variance of every variable to be 1. 
# this should not be all that difficult given that the confounders are exogenous variables, however,
# we need to make sure that we do not 'blow' up the effect of the confounders to such an extend that there is no longer
# an innovation matrix that makes the system stationary. 

set.seed(427)
library(mvtnorm)

N <- 100000
T <- 5

ax <- 0.25
ay <- 0.25
bx <- 0.10
by <- 0.10
rho <- 0.30

A <- matrix(c(ax, by, bx, ay), 2, 2, byrow=TRUE)

k <- 3
tau2 <- rep(1, k)
Psi <- diag(tau2)

# time-varying B matrices
# meaning we specify a list rather than a single matrix
Beta_list <- list(
  rbind(c(0.30,0.20,0.10), c(0.30,0.20,0.10)),
  rbind(c(0.25,0.18,0.08), c(0.25,0.18,0.08)),
  rbind(c(0.20,0.16,0.06), c(0.20,0.16,0.06)),
  rbind(c(0.15,0.14,0.04), c(0.15,0.14,0.04)),
  rbind(c(0.10,0.12,0.02), c(0.10,0.12,0.02))
)

find_c <- function(A, S_U, rho) {
  f <- function(c) {
    S_target_c <- matrix(c(1,c,c,1),2,2)
    S_dyn_c <- S_target_c - S_U
    Sigma_e_c <- S_dyn_c - t(A) %*% S_dyn_c %*% A
    v1 <- Sigma_e_c[1,1]
    v2 <- Sigma_e_c[2,2]
    cov12 <- Sigma_e_c[1,2]
    cov12/sqrt(v1*v2) - rho
  }
  uniroot(f, c(-0.99,0.99))$root
}

# compute mean confounder-variance contribution across waves
S_U_list <- lapply(Beta_list, function(B) B %*% Psi %*% t(B))
S_U <- Reduce("+", S_U_list) / length(S_U_list)

c_stat <- find_c(A, S_U, rho)

S_target <- matrix(c(1,c_stat,c_stat,1),2,2)
S_dyn <- S_target - S_U
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2
Sigma_e1 <- (S_dyn + t(S_dyn))/2

U <- rmvnorm(N, mean=rep(0,k), sigma=Psi)

df <- matrix(NA, nrow=N, ncol=2*T + k)
colnames(df) <- c(paste0("x",1:T), paste0("y",1:T), paste0("c",1:k))
df[, (2*T+1):(2*T+k)] <- U

Ddyn <- rmvnorm(N, mean=c(0,0), sigma=Sigma_e1)
obs1 <- Ddyn + U %*% t(Beta_list[[1]])
df[,1] <- obs1[,1]
df[,1+T] <- obs1[,2]

for(i in 2:T){
  Ddyn <- Ddyn %*% t(A) + rmvnorm(N, sigma=Sigma_e)
  # use B_t for each wave
  obs <- Ddyn + U %*% t(Beta_list[[i]])
  df[,i] <- obs[,1]
  df[,i+T] <- obs[,2]
}

df <- as.data.frame(df)



