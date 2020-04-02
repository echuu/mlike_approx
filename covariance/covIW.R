
## covIW.R

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("covariance/covIW_helper.R")

library(matrixcalc) # move into one of the other header files later
library(mvtnorm)

## define necessary parameters needed for this simulation

J = 500                # num of MCMC samples (from true posterior in this case)
N = 500                # number of observations
p = 5                  # num rows/cols in the covariance matrix
D = 0.5 * p * (p + 1)  # dimension of u that is fed into the tree


## wishart prior parameters
Omega = diag(1, p)  # scale matrix
nu    = p + 1       # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(rWishart(1, p, Omega), p)
is.positive.definite(Sigma)

## generate data
X = rmvnorm(N, mean = rep(0, p), sigma = Sigma) # (N x p)
S = matrix(0, p, p)
for (n in 1:N) {
    S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
}

# ------------------------------------------------------------------------------

## TODO: obtain posterior samples
Sigma_post = solve(Wishart_InvA_RNG(nu + N, S + Omega))
Sigma_post
Sigma


## TODO: fill out the helper functions
## TODO: modify algorithm so that only constant approximations are made
## TODO: compare approximations to true log marginal likelihood
## TODO: perform asymptotic analysis
















