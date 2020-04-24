
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
D = 5                       # num rows/cols in the covariance matrix
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


## (3a) compute maximized likelihood (using true Sigma)
loglik_max = maxLogLik(Sigma, param_list)


## (3b) compute approximation
hml_approx = hml_const(1, D_u, u_df, J, param_list)


# the log ML approximation is stored in the "const_vec" variable
# subtract off the maximized log likelihood
hml_approx$const_vec - loglik_max


# (3c) compute true log ML, subtract off maximized log likelihood
lil(param_list) - maxLogLik(Sigma, param_list)

hml_approx$const_vec
lil(param_list)

# ------------------------------------------------------------------------------

#### model output diagnostics

part_info = hml_approx$param_out
u_df_info = hml_approx$u_df_fit

# (1.1) for each point: fitted value, leaf_id, residual

## notation: 
# psi_u     : true value, psi(u)
# psi_star  : psi(u_star), where u_star is the partition's representative point
# psi_resid : psi_u - psi_star
# psi_hat   : fitted value for psi(u) on a given partition
##

u_df_info %>% head

# (1.2) for each partition: 'median' points, fitted value, upper/lower bounds,
# number of observations in partition
part_info 


# compare the psi values from tree vs. psi values evaluated at 'representative'
# point of each partition
part_psi = part_info %>% 
    dplyr::select(leaf_id, psi_hat) %>% 
    merge(u_df_info %>% dplyr::select(leaf_id, psi_star) %>% unique, 
          by = 'leaf_id')

part_psi

# (2) fitted value for each partition (const_approx)
# (2.1) look at the MSE for each of the partitions to see if there's one
# parition where there is significantly larger discrepancy
# (2.2) consider the number of observations in each of the partitions in
# conjunction with the MSE for each partition

library(MLmetrics)
part_mse = u_df_info %>% 
    dplyr::group_by(leaf_id) %>% 
    summarise(psi_star_mse = MSE(psi_u, psi_star)) %>% 
    merge(part_info %>% dplyr::select(leaf_id, n_obs), by = 'leaf_id')


# for each partition, display the fitted value for psi (from tree), 
# psi evaluated at the representative point, the mse associated with the
# representative point, the number of observations for that partition
psi_df = merge(part_psi, part_mse, by = 'leaf_id')

psi_df %>% arrange(psi_star_mse)

u_df_info %>% filter(leaf_id == 28)

psi_df %>% arrange(n_obs)

# compute mse for the tree's fitted values
psi_rpart = hml_approx$u_rpart
tree_dev = rpart_mse(psi_rpart, u_df)
tree_mse = tree_dev %>% group_by(leaf_id) %>% 
    summarise(psi_hat_mse = mean(dev_sq))

# gather all psi approximations (tree and algo) with associated MSEs
# into one table
psi_mse = psi_df %>% merge(tree_mse, by = 'leaf_id') %>% 
    dplyr::select(leaf_id, psi_hat, psi_hat_mse, psi_star, psi_star_mse, n_obs)


psi_mse %>% arrange(n_obs) %>% mutate(n_perc = n_obs / J)


hml_approx$const_vec
lil(param_list)



# (3) compute approximate integral over each partition
hml_approx$const_approx


# (4) extract partitions corresponding to the diagonal

# extract column names corresponding to diagonal entries (this is built into 
# getDiagCols() function already, but could be of some use later to have these
# indices available
diagInd = getDiagIndex(D, D_u)
diag_names = paste("u", diagInd, sep = '')

# extract u_k_star with corresponding lb, ub for each k that is a diagonal entry
getDiagCols(part_info, D, D_u)













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











