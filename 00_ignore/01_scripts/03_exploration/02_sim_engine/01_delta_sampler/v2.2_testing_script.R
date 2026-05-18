# goal is to stress test all parts of the data generation process
# to make sure that everything is working as expected.
# also we have a look at the helper functions to understand how they work. 

###########################################################################
# DELTA SAMPLER
###########################################################################

# first we start by looking at the helper functions
# function to make the index names for two-way interactions
# --------------------------------------------------------------------------
all_pairs <- function(k) {
  t(combn(seq_len(k), 2))
}

# test that the function works as expected
all_pairs(5)

# function to make the index names for three-way interactions.
# --------------------------------------------------------------------------
all_triples <- function(k) {
  t(combn(seq_len(k), 3))
}

# test that the function works as expected
all_triples(5)

# function to compute all possible pairings of 6 elements
# we need those to finally compute E[C1 C2 C3 C4 C5 C6]
# this could also be done by storing all combinations, as this is a fixed set 
# --------------------------------------------------------------------------
get_pairings6 <- function() {
  rec <- function(v) {
    if (length(v) == 0) return(list(list()))
    if (length(v) == 2) return(list(list(c(v[1], v[2]))))

    first <- v[1]
    out <- list()

    for (m in 2:length(v)) {
      rest <- v[-c(1, m)]
      sub  <- rec(rest)

      for (s in sub) {
        out[[length(out) + 1]] <- c(list(c(first, v[m])), s)
      }
    }

    out
  }
  rec(1:6)
}

# test that the function works as expected
get_pairings6()

# then the function to actually compute E[C1 C2 C3 C4 C5 C6]
# --------------------------------------------------------------------------
sixth_moment_gaussian <- function(idx, Omega11, pairings6 = get_pairings6()) {
  sum(vapply(
    pairings6,
    function(pr) {
      prod(vapply(
        pr,
        function(pa) Omega11[idx[pa[1]], idx[pa[2]]],
        numeric(1)
      ))
    },
    numeric(1)
  ))
}

# test that the function works as expected
pairings6 <- get_pairings6()

Omega11 <- matrix(c(
  1.0, 0.3, 0.2,
  0.3, 1.0, 0.4,
  0.2, 0.4, 1.0
), nrow = 3, byrow = TRUE)

Omega11_1 <- diag(3)

sixth_moment_gaussian(c(1, 2, 3, 1, 2, 3), Omega11, pairings6)

# now to check the result:
# extract correlations
rho12 <- Omega11[1,2]
rho13 <- Omega11[1,3]
rho23 <- Omega11[2,3]

# closed-form formula
  1 +
  2*rho12^2 +
  2*rho13^2 +
  2*rho23^2 +
  8*rho12*rho13*rho23

# function to compute the covariance between two standardized two-way interactions
# --------------------------------------------------------------------------
cov_Z2 <- function(i, j, p, q, Omega11) {
  num <- Omega11[i, p] * Omega11[j, q] +
         Omega11[i, q] * Omega11[j, p]

  den <- sqrt((1 + Omega11[i, j]^2) *
              (1 + Omega11[p, q]^2))

  num / den
}

# Example:
# compute Cov(Z_12, Z_13)
i <- 1
j <- 2
p <- 1
q <- 3

cov_Z2(i, j, p, q, Omega11)
cov_Z2(i, j, p, q, Omega11_1)

# function to compute the covariance between two standardized three-way interactions
# --------------------------------------------------------------------------
cov_T3 <- function(i, j, l, p, q, r, Omega11, pairings6 = get_pairings6()) {
  num <- sixth_moment_gaussian(c(i, j, l, p, q, r), Omega11, pairings6)

  den1 <- 1 +
    2 * Omega11[i, j]^2 +
    2 * Omega11[i, l]^2 +
    2 * Omega11[j, l]^2 +
    8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]

  den2 <- 1 +
    2 * Omega11[p, q]^2 +
    2 * Omega11[p, r]^2 +
    2 * Omega11[q, r]^2 +
    8 * Omega11[p, q] * Omega11[p, r] * Omega11[q, r]

  num / sqrt(den1 * den2)
}

# example
i <- 1
j <- 2
l <- 3
p <- 1
q <- 2
r <- 3

# the covariance between T_123 and T_123 = var(T_123)
cov_T3(i, j, l, p, q, r, Omega11)
cov_T3(i, j, l, p, q, r, Omega11_1)

# covariance between T_123 and T_234
Omega11_4 <- matrix(c(
  1.0, 0.3, 0.2, 0.1,
  0.3, 1.0, 0.4, 0.2,
  0.2, 0.4, 1.0, 0.3,
  0.1, 0.2, 0.3, 1.0
), nrow = 4, byrow = TRUE)

i <- 1
j <- 2
l <- 3
p <- 2
q <- 3
r <- 4

cov_T3(i, j, l, p, q, r, Omega11_4)

# function to compute the covariance between main effect and three-way interaction
# --------------------------------------------------------------------------
cov_13 <- function(a, i, j, l, Omega11) {
  num <- Omega11[a, i] * Omega11[j, l] +
         Omega11[a, j] * Omega11[i, l] +
         Omega11[a, l] * Omega11[i, j]

  den <- sqrt(
    1 +
    2 * Omega11[i, j]^2 +
    2 * Omega11[i, l]^2 +
    2 * Omega11[j, l]^2 +
    8 * Omega11[i, j] * Omega11[i, l] * Omega11[j, l]
  )

  num / den
}
# test Cov(C1, T234)
a <- 1
i <- 2
j <- 3
l <- 4

cov_13(a, i, j, l, Omega11_4)

l <- 1
cov_13(a, i, j, l, Omega11_1)
cov_13(a, i, j, l, Omega11)

# function to build the varcov matrix of two-way interactions
# --------------------------------------------------------------------------
build_Omega22 <- function(Omega11) {
  P <- all_pairs(nrow(Omega11))
  m2 <- nrow(P)

  Omega22 <- matrix(0, m2, m2)
  if (m2 == 0) return(Omega22)

  diag(Omega22) <- 1

  for (a in seq_len(m2)) {
    i <- P[a, 1]
    j <- P[a, 2]

    for (b in seq_len(m2)) {
      if (a == b) next

      p <- P[b, 1]
      q <- P[b, 2]

      Omega22[a, b] <- cov_Z2(i, j, p, q, Omega11)
    }
  }

  Omega22
}

build_Omega22(Omega11)
build_Omega22(Omega11_1)

# function to build the three-way interaction covariance matrix
# --------------------------------------------------------------------------
build_Omega33 <- function(Omega11, pairings6 = get_pairings6()) {
  T <- all_triples(nrow(Omega11))
  m3 <- nrow(T)

  Omega33 <- matrix(0, m3, m3)
  if (m3 == 0) return(Omega33)

  diag(Omega33) <- 1

  for (a in seq_len(m3)) {
    i <- T[a, 1]
    j <- T[a, 2]
    l <- T[a, 3]

    for (b in seq_len(m3)) {
      if (a == b) next

      p <- T[b, 1]
      q <- T[b, 2]
      r <- T[b, 3]

      Omega33[a, b] <- cov_T3(i, j, l, p, q, r, Omega11, pairings6)
    }
  }

  Omega33
}

# example
build_Omega33(Omega11_4)
build_Omega33(Omega11_1)

# function to build the main with three-way interaction covariance matrix
# --------------------------------------------------------------------------
build_Omega13 <- function(Omega11, k) {
  T <- all_triples(nrow(Omega11))
  m3 <- nrow(T)

  Omega13 <- matrix(0, k, m3)
  if (m3 == 0) return(Omega13)

  for (t in seq_len(m3)) {
    i <- T[t, 1]
    j <- T[t, 2]
    l <- T[t, 3]

    for (a in seq_len(k)) {
      Omega13[a, t] <- cov_13(a, i, j, l, Omega11)
    }
  }

  Omega13
}

# example
build_Omega13(Omega11_4, 4)
build_Omega13(Omega11_1, 3)

# check the functionality of the entire function
library(here)
source(here("01_scripts", "03_exploration", "02_sim_engine", "01_delta_sampler", "v2_delta_sampler.R"))

out_1 <- sample_delta_1(
  k = 1,
  Omega11 = diag(1),
  R2_total = 0.20,
  R2_interaction = 0,
  include_2way = FALSE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

out_2 <- sample_delta_1(
  k = 2,
  Omega11 = diag(2),
  R2_total = 0.20,
  R2_interaction = 0.1,
  include_2way = TRUE,
  include_3way = FALSE,
  min_abs = 0,
  max_abs = 1
)

out_2$Delta

out_3 <- sample_delta_1(
  k=3,
  Omega11 = Omega11_1,
  R2_total = 0.15,
  R2_interaction = 0.3,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

out_3$Delta

out_4 <- sample_delta_1(
  k=4,
  Omega11 = Omega11_4,
  R2_total = 0.15,
  R2_interaction = 0.3,
  include_2way = TRUE,
  include_3way = TRUE,
  min_abs = 0,
  max_abs = 1
)

out_4$Delta
