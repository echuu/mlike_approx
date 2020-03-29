

library(dplyr)
library(mvtnorm)

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx_v1.R")            # load main algorithm functions
source("RRR/RRR_C_helper.R")            # load psi(), lambda(), preprocess()
source("RRR/gibbs_RRR.R")               # load sampleRRR()
source('extractPartition.R')            # load extractPartition() function


p = 3             # number of columns in X
q = 5             # number of columns in Y

r_0 = 2           # true rank
r = 2             # model-postulated rank
D = r * p + q * r # dimension of each MCMC sample
D_C = p * q
sig2 = 1          # fixed for now.
del = 10^(-2)     # prior parameter

set.seed(1)
A_0 = matrix(rnorm(p * r_0, 0, 1), p, r_0) # (p x r) matrix
B_0 = matrix(rnorm(q * r_0, 0, 1), q, r_0) # (q x r) matrix
C_0 = A_0 %*% t(B_0)

nMCMC = 300       # number of MCMC samples from the posterior AFTER burnin
nBurn = 500       # number of samples to discard


dim_model = r_0 * (p + q - r_0)
    
(true_slope = dim_model / 2) # 6 for p = 3, q = 5, r_0 = 2, r = 2

# test one case ----------------------------------------------------------------

N = 100

gibbs_obj = sampleRRR(nMCMC, nBurn, A_0, B_0, p, q, r, r, D, N, sig2, del)

param_list = list(p = p, q = q, r = r, n = N, d = D,          # dimensions vars
                  Y = gibbs_obj$Y, X = gibbs_obj$X,           # response, design
                  XtX = gibbs_obj$XtX, Xty = gibbs_obj$Xty,   # precompute
                  sig2 = sig2, del = del)                     # prior params

u_samps = gibbs_obj$u_samps
c_samps = gibbs_obj$c_samps

c_samps[300,] %>% matrix(p, q) # this should be 'close' to C_0

u = u_samps[300,] %>% unname %>% unlist
A_post = u[1:(p*r)] %>% matrix(p, r)
Bt_post = u[(p * r + 1):D] %>% matrix(r, q)

A_post %*% Bt_post             # this should equal the line below
c_samps[300,] %>% matrix(p, q) # this should be 'close' to C_0

c_df = preprocess(c_samps, D_C, param_list)

hml_approx = hml(1, D_C, c_df, nMCMC, param_list)
hml_approx$hybrid_vec

# ------------------------------------------------------------------------------

K_sims    = 100                                # num of sims to run for each N
N_vec_log = seq(5, 9, by = 0.2)              # sample size grid unif in log
N_vec     = floor(exp(N_vec_log)) %>% unique   # sample size to generate data
LIL_N_k_hat = matrix(0, length(N_vec), K_sims) 


for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    for (k in 1:K_sims) {
        
        # generate data, sample from posterior ---------------------------------
        gibbs_obj = sampleRRR(nMCMC, nBurn, A_0, B_0, p, q, r_0, r, D, N, sig2, del)
        
        param_list = list(p = p, q = q, r = r, n = N, d = D,  # dimensions vars
                          Y = gibbs_obj$Y, X = gibbs_obj$X,   # response, design
                          XtX = gibbs_obj$XtX, Xty = gibbs_obj$Xty,
                          sig2 = sig2, del = del)             # prior params
        
        # compute log-likelihood evaluated at max
        loglik_max = loglik_true(A_0, B_0, param_list)
        
        # extract posterior samples from gibbs object
        c_samps = gibbs_obj$c_samps
        
        # evaluate psi(u) for each of the posterior samples
        c_df = preprocess(c_samps, D_C, param_list) # nMCMC x (d + 1) 
        
        # generate hybrid approximation
        hml_approx = hml(1, D_C, c_df, nMCMC, param_list)
        
        # subtract maximized likelihood from the resulting approximation
        LIL_N_k_hat[i, k] = hml_approx$hybrid_vec - loglik_max # -147.6771
        
    } # end inner loop over K_sims
    
    print(paste("iter ", i, "/", length(N_vec), ": ",
                "approx LIL for N = ", N, " -- LIL = ",
                round(mean(LIL_N_k_hat[i, ]), 2), sep = ''))
    
} # end outer loop over sample size
















