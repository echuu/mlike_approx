
# important: files below must be sourced in order
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("covarIW_helper.R")  # covariance related helper functions

library(Rcpp)

sourceCpp("C:/Users/ericc/mlike_approx/speedup/fast_covIW.cpp")


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

# log_density = function(u, data) {
#     -psi(u, data)
# }


# ## (2) obtain posterior samples
# postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post
# 
# # these are the posterior samples stored as vectors (the lower cholesky factors
# # have been collapsed into D_u dimensional vectors)
# post_samps = postIW$post_samps                 # (J x D_u)
# u_samp = as.matrix(post_samps)
# colnames(u_samp) = names(u_df)[1:D_u]
# # prepare bridge_sampler input()
# lb = rep(-Inf, D_u)
# ub = rep(Inf, D_u)
# 
# diag_ind = getDiagIndex(D, D_u) # obtain column index of the diagonal entries
# lb[diag_ind] = 0                # diagonal entries are positive
# names(lb) <- names(ub) <- colnames(u_samp)
# 
# bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
#                                                log_posterior = log_density,
#                                                data = param_list, 
#                                                lb = lb, ub = ub, 
#                                                silent = TRUE)
# bridge_result$logml

(true_logml = lil(param_list)) 

B = 10 # number of replications
hyb_wt1  = numeric(B)  # store harmonic mean estiator
hyb_wt2  = numeric(B)  # store harmonic mean estiator
hyb_avg  = numeric(B)  # store harmonic mean estiator
# bridge = numeric(B)
set.seed(1)
for (b_i in 1:B) {
    
    postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post
    
    # these are the posterior samples stored as vectors (the lower cholesky factors
    # have been collapsed into D_u dimensional vectors)
    post_samps = postIW$post_samps                 # (J x D_u)
    
    
    # u_df stores the posterior samples row-wise so that the first D_u columns 
    # store the lower cholesky factors in vector form, and the last column is
    # the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
    u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
    #
    # ## (3b) compute approximation
    hybrid = logml(D_u, u_df, J, param_list)
    
    if (any(is.na(hybrid))) {print(paste("error in iteration", b_i)); next;}
    
    hybrid$all_approx
    hybrid$wt_approx1
    hybrid$wt_approx2
    hybrid$wt_approx3
    
    hyb_wt1[b_i] = hybrid$wt_approx1
    hyb_wt2[b_i] = hybrid$wt_approx2
    hyb_wt3[b_i] = hybrid$wt_approx3

    # u_samp = as.matrix(post_samps)
    # colnames(u_samp) = names(u_df)[1:D_u]
    # 
    # bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
    #                                                log_posterior = log_density,
    #                                                data = param_list, 
    #                                                lb = lb, ub = ub, 
    #                                                silent = TRUE)
    # bridge[b_i] = bridge_result$logml
    # 
    # if (b_i %% 10 == 0) { 
    #     print(paste("iter: ", b_i, ' - ',
    #                 round(mean(bridge[1:b_i]), 3), 
    #                 sep = '')) 
    # }
    
}

LIL = true_logml
approx = data.frame(bridge)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)









