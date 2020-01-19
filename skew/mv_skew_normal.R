



library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function
library('microbenchmark')

# path for lenovo
setwd("C:/Users/ericc/mlike_approx")

# path for dell
# setwd("C:/Users/chuu/mlike_approx")
source("partition/partition.R")      # load partition extraction functions

library(sn)
library(VGAM)

D = 2
N = 50 # pseudo-sample size

set.seed(1)

J = 200

Omega = diag(1, D)
Sigma = D / N * Omega   # to mimic situation where gamma is relatively concentrated
Sigma_inv = solve(Sigma)
alpha = c(1, 2)  # slant

mu_0 = c(5, 2)

u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame # 200 x 3

# ------------------------------------------------------------------------------


u_df_full = preprocess(u_samps, D)

approx_skew = approx_lil_stan(1, D, u_df_full, J)

approx_skew




D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) # -1.4633 for D = 3


# ------------------------------------------------------------------------------
D = 2
N = 50 # pseudo-sample size

set.seed(1)

J = 2000

Omega = diag(1, D)
Sigma = D / N * Omega   # to mimic situation where gamma is relatively concentrated
Sigma_inv = solve(Sigma)
alpha = c(1, 2)  # slant

mu_0 = c(0, 0)

u_samps = rmvnorm(J, mu_0, sigma = Sigma) %>% data.frame # 200 x 3


# ------------------------------------------------------------------------------


u_df_full = preprocess(u_samps, D)

approx_skew = approx_lil_stan(1, D, u_df_full, J)

log(-sum(approx_skew))

D / 2 * log(2 * pi) + 0.5 * log_det(Sigma)


D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) # -1.4633 for D = 3







