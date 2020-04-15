
# covarIW.R 


# change to your working directory
# note: all files must be in the same directory
# setwd("C:/Users/ericc/Dropbox/logML")

# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("covarIW_helper.R")  # covariance related helper functions


N = 100                     # number of observations
D = 3                       # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree
J = 300

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
S = matrix(0, D, D)
for (n in 1:N) {
    S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
}


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


## (3a) compute maximized likelihood (using true Sigma)
loglik_max = maxLogLik(Sigma, param_list)


## (3b) compute approximation
hml_approx = hml_const(1, D_u, u_df, J, param_list)


# the log ML approximation is stored in the "const_vec" variable
# subtract off the maximized log likelihood
hml_approx$const_vec - loglik_max


# (3c) compute true log ML, subtract off maximized log likelihood
lil(param_list) - maxLogLik(Sigma, param_list)



# ------------------------------------------------------------------------------

library(cubature)
library(matrixcalc)


e_psi = function(u) {
    
    exp(-psi(u, param_list))
    
}

u_0 = u_df[1,-4] %>% unname %>% unlist
L_test = matrix(0, D, D)
L_test[lower.tri(L_test, diag = T)] = u_0

stable_numer = adaptIntegrate(e_psi, 
                              lowerLimit = c(0, -Inf, -Inf, 0, -Inf, 0), 
                              upperLimit = rep(Inf, D_u),
                              tol = 1e-4)

hml_approx$const_vec
log(stable_numer$integral)
lil(param_list)











