# this script solves the problem from the previous script by simulating effects based on some random walk
# we need to think about ways that the system does not blow up
# one of those ways is to have some shrinkage toward a baseline effect, which will act as soft bounds

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
Psi <- diag(rep(1,k))

B0 <- rbind(c(0.30,0.20,0.10),
            c(0.30,0.20,0.10))

# RW mean vector for each coefficient (same dimension as B0)
# basically denoting overtime trends in the effects 
RW_mean <- rbind(c(0.00, 0.00, 0.00),
                 c(0.00, 0.00, 0.00))

# now we set the random walk sd and shrinkage toward B0
# the sd denotes how much variability there is in the random walk step
# the shrinkage denotes how much pull there is back toward B0 at each time step
rw_sd  <- 0.03
shrink <- 0.10

Beta_list <- vector("list", T)
Beta_list[[1]] <- B0

# function to create random walk with shrinkage toward B0
for(t in 2:T){

  # random walk step with mean vector
  step <- matrix(
    rnorm(length(B0), mean = as.vector(RW_mean), sd = rw_sd),
    nrow = 2
  )

  B_new <- Beta_list[[t-1]] + step + shrink*(B0 - Beta_list[[t-1]])

  Beta_list[[t]] <- B_new
}

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

S_U_list <- lapply(Beta_list, function(B) B %*% Psi %*% t(B))
S_U <- Reduce("+", S_U_list) / length(S_U_list)

c_stat <- find_c(A, S_U, rho)

S_target <- matrix(c(1,c_stat,c_stat,1),2,2)
S_dyn <- S_target - S_U
Sigma_e <- S_dyn - t(A) %*% S_dyn %*% A
Sigma_e <- (Sigma_e + t(Sigma_e))/2
Sigma_e1 <- S_dyn

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
  obs <- Ddyn + U %*% t(Beta_list[[i]])
  df[,i] <- obs[,1]
  df[,i+T] <- obs[,2]
}

df <- as.data.frame(df)
df
