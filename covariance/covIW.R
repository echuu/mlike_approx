
## covIW.R

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx_v1.R")
source("covariance/covIW_helper.R")

library(matrixcalc) # move into one of the other header files later
library(mvtnorm)

## define necessary parameters needed for this simulation

J = 500                  # num of MCMC samples (from true posterior in this case)
N = 500                  # number of observations
D = 5                    # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)  # dimension of u that is fed into the tree


## wishart prior parameters
Omega = diag(1, D)  # scale matrix
nu    = D + 1       # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(rWishart(1, D, Omega), D)
is.positive.definite(Sigma)

## generate data
X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = matrix(0, D, D)
for (n in 1:N) {
    S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
}

# ------------------------------------------------------------------------------

param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu,         # prior parameters
                  u_df = NULL)                    # posterior samples

## obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega) # post_samps, Sigma_post, L_post

Sigma0 = postIW$Sigma_post[[1]]

u = postIW$post_samps[1,]
L = postIW$L_post[[1]]

L %*% t(L)
Sigma0

cov_logprior(u, param_list)             # -47.20618
cov_logprior_sigma(Sigma0, param_list)
cov_logprior_L(L, param_list)

L1 = matrix(0, D, D)
L1[lower.tri(L1, diag = T)] = postIW$post_samps[1,]






postIW$Sigma_post[[1]] %>% det


(postIW$L_post[[1]] %>% det)^2

(postIW$L_post[[1]] %>% diag %>% prod)^2

## TODO: fill out the helper functions
## TODO: modify algorithm so that only constant approximations are made
## TODO: compare approximations to true log marginal likelihood
## TODO: perform asymptotic analysis



## verify that we have the correct expression for the true logML by plotting
## logML vs. logN for 50 replicates per N

J           = 300                               # num of MCMC samples from post
K_sims      = 50                                # num of sims to run for each N
N_vec_log   = seq(5, 12, by = 0.25)             # sample size grid unif in log
N_vec       = floor(exp(N_vec_log)) %>% unique  # sample size to generate data
LIL_N_k_hat = matrix(0, length(N_vec), K_sims)  # store approximations
LIL_N       = numeric(length(N_vec))            # store true logML

length(N_vec)


for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    LIL_N_k = numeric(K_sims)
    
    for (k in 1:K_sims) {
        
        X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
        S = matrix(0, D, D)
        for (n in 1:N) {
            S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
        }
        
        param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                          Omega = Omega, nu = nu,         # prior parameters
                          u_df = NULL)                    # posterior samples
        LIL_N_k[k] = lil(param_list) - maxLogLik(Sigma, param_list)
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
}


LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))
# write.csv(LIL_df, "true_lil.csv", row.names = F) # N = seq(5, 12, by = 0.25)


library(reshape2)
library(ggpmisc)

formula1 = y ~ x


ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Hybrid (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8)


# ------------------------------------------------------------------------------

## testing collapse functionality
D = 3
Sigma = matrix(rWishart(1, D, diag(1, D)), D)

L = t(chol(Sigma))

u = L[lower.tri(L, diag = T)]     # map lower diagonal + diagonal to a vector 

L1 = matrix(0, D, D)
L1[lower.tri(L1, diag = T)] = u   # map vector back to lower diagonal matrix














