
## covIW.R



source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx_v1.R")
source("hybrid_const.R")
source("covariance/covIW_helper.R")

library(matrixcalc) # move into one of the other header files later
library(mvtnorm)

library(MCMCpack)
library(CholWishart)


options(scipen = 999) # moved to load_libraries.R file


# ------------------------------------------------------------------------------



setwd("C:/Users/ericc/Dropbox/logML") # change to your working directory
source("setup.R")

LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)
source("covariance/covIW_helper.R")


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

# ------------------------------------------------------------------------------

## generate data
X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = matrix(0, D, D)
for (n in 1:N) {
    S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
}

# ------------------------------------------------------------------------------

param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                  Omega = Omega, nu = nu)         # prior parmas

## obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega) # post_samps, Sigma_post, L_post


Sigma_post = riwish(nu + N, Omega + S)
Sigma_post

(Sigma0 = postIW$Sigma_post[[1]])
Sigma


L = t(chol(Sigma))

L %*% t(L)
Sigma0


L1 = matrix(0, D, D)
L1[lower.tri(L1, diag = T)] = postIW$post_samps[1,]


# cov_logprior_sigma(Sigma0, param_list)
# cov_logprior_L(L, param_list)

# test log-prior function
cov_logprior(u, param_list)       # -34.36231

# test log-likelihood function
cov_loglik(u, param_list)

u = postIW$post_samps[5,]
# test psi() function
psi(u, param_list)

# evaluate psi() for each of the rows in post_samps
post_samps = postIW$post_samps                   # (J x D_u)
u_df = preprocess(post_samps, D_u, param_list)   # J x (D_u + 1)



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

i = 1; k = 1;

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    
    LIL_N_k = numeric(K_sims) # store true log ML
    
    for (k in 1:K_sims) {
        set.seed(1)
        X = rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
        S = matrix(0, D, D)
        for (n in 1:N) {
            S = S + tcrossprod(X[n,]) # compute sum_n x_n * x_n'
        }
        
        param_list = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
                          Omega = Omega, nu = nu,         # prior parameters
                          u_df = NULL)                    # posterior samples
        
        # compute maximized log likelihood
        loglik_max = maxLogLik(Sigma, param_list)
        
        # postIW contains: post_samps, Sigma_post, L_post
        postIW = sampleIW(J, N, D_u, nu, S, Omega) 
        
        post_samps = postIW$post_samps                   # (J x D_u)
        u_df = preprocess(post_samps, D_u, param_list)   # J x (D_u + 1)
    
        # generate hybrid approximation
        # hml_approx() is a version of hml() that ignores the gradient term
        hml_approx = hml_const(1, D_u, u_df, J, param_list)
        
        # subtract maximized likelihood from the resulting approximation
        LIL_N_k_hat[i, k] = hml_approx$const_vec - loglik_max
        
        hml_approx$const_vec - loglik_max # -34.70373
        
        # compute true log ML - maximized likelihood
        LIL_N_k[k] = lil(param_list) - maxLogLik(Sigma, param_list)
        
        lil(param_list) - maxLogLik(Sigma, param_list) # -43.2191
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
    print(paste("iter ", i, "/", length(N_vec), ": ",
                "approx LIL for N = ", N, " -- LIL = ",
                round(mean(LIL_N_k_hat[i, ]), 2), 
                " (", round(LIL_N[i], 2), ")", 
                sep = ''))
    
}


lil_hyb   = rowMeans(LIL_N_k_hat)  # length(N_vec) x 1
log_N     = log(N_vec)             # length(N_vec) x 1

LIL_df = data.frame(LIL_hat = lil_hyb, log_N = log_N)
# LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))
# write.csv(LIL_df, "true_lil.csv", row.names = F) # N = seq(5, 12, by = 0.25)


LIL_df = read.csv("covariance/true_lil.csv")

library(reshape2)
library(ggpmisc)

formula1 = y ~ x


# approx
ggplot(LIL_df, aes(x = log_N, y = LIL_hat)) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Hybrid (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8)


# true
ggplot(LIL_df, aes(x = log_N, y = LIL_N)) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Hybrid (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8)


# TODO: with the true log marginal likelihood
library(reshape2)
lil_0   = LIL_N                  # length(N_vec) x 1
lil_hyb = rowMeans(LIL_N_k_hat)  # length(N_vec) x 1
log_N   = log(N_vec)             # length(N_vec) x 1

LIL_df = data.frame(LIL_N = lil_0, LIL_hat = lil_hyb, log_N = log(N_vec))

LIL_df_long = melt(LIL_df, id.vars = "log_N")
head(LIL_df_long)

ggplot(LIL_df_long, aes(x = log_N, y = value, 
                        color = as.factor(variable))) + geom_point(size = 1.3) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "True (Red), Approx (Blue)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")



# ------------------------------------------------------------------------------

## testing collapse functionality
D = 3
Sigma = matrix(rWishart(1, D, diag(1, D)), D)

L = t(chol(Sigma))

u = L[lower.tri(L, diag = T)]     # map lower diagonal + diagonal to a vector 

L1 = matrix(0, D, D)
L1[lower.tri(L1, diag = T)] = u   # map vector back to lower diagonal matrix














