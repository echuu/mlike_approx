

# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")

## look in here for psi implementation
source("covarIW_helper.R")  # covariance related helper functions


N = 100                     # number of observations
D = 8                       # num rows/cols in the covariance matrix


## in this example D_u plays the role that D played in previous examples
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree


J = 1000


## wishart prior parameters
Omega = diag(1, D)          # scale matrix
nu    = D + 1               # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(rWishart(1, D, Omega), D)
is.positive.definite(Sigma)


##
## below is one iteration of the simulation: 
##     (1) generate data
##     (2) sample from posterior
##     (3) compute: 
##                  (a) maximized likelihood 
##                  (b) approximate logML
##                  (c) true logML
## 

## (1) generate data
X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = t(X) %*% X                                  # (p x p)


## store parameters in a list that can be passed into the algorithm
param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu)         # prior params


## (2) obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post


# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)


# u_df stores the posterior samples row-wise so that the first D_u columns 
# store the lower cholesky factors in vector form, and the last column is
# the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)

head(u_df)


hml_approx = hml_const(1, D_u, u_df, J, param_list)


hml_approx$const_vec # check this after doing fast/slow psi
