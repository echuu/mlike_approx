

# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("covarIW_helper.R")  # covariance related helper functions

N = 100                     # number of observations
D = 5                       # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree
J = 5000


## wishart prior parameters
Omega = diag(1, D)          # scale matrix
nu    = D + 1               # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(rWishart(1, D, Omega), D)
is.positive.definite(Sigma)


## (1) generate data
X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = t(X) %*% X                                  # (p x p)


## store parameters in a list that can be passed into the algorithm
param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu)         # prior params


postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post

# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)



## slow version of psi
u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)


## fast version of psi
sourceCpp("C:/Users/ericc/mlike_approx/speedup/fast_covIW.cpp")
u_df_fast = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)


# compare psi values in the following: 
u_df$psi_u %>% head

u_df_fast$psi_u %>% head






