

# ------------------------------------------------------------------------------


library(mvtnorm)           # for draws from multivariate normal
library('MCMCpack')        # for rinvgamma() function
library(sn)
library(VGAM)


DELL_PATH = "C:/Users/chuu/mlike_approx"
# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
# setwd(LEN_PATH)

# path for dell
setwd(DELL_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions

source("skew/mv_skew_normal_helper.R")  # load psi(), lambda()


# fixed settings ---------------------------------------------------------------
D = 2
N = 500 # pseudo-sample size
Omega = diag(1, D)
Sigma = D / N * Omega 
Sigma_inv = solve(Sigma)
alpha = rep(1, D) 
mu_0 = rep(0, D)

prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, 
             alpha = alpha, mu_0 = mu_0)

# -4.3767 for D = 2, N = 500
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) 

# ------------------------------------------------------------------------------

set.seed(1)
J = 5000
N_approx = 10
u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J/N_approx, prior)
mean(approx_skew) # -3.085342


# ------------------------------------------------------------------------------
plot(u_samps)
u_df_full = preprocess(u_samps, D)

# skew_rpart = rpart(psi_u ~ ., u_df_full)
# plot(skew_rpart)
# approx_skew

D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) # -4.376731 for D = 2

# ------------------------------------------------------------------------------


# repeat analysis for many D
set.seed(1)
N = 5000 # pseudo-sample size


D_vec = c(2:10)
LIL_D = numeric(length(D_vec))
LIL_D_approx = numeric(length(D_vec))
for (d_i in 1:length(D_vec)) {
    
    D = D_vec[d_i]
    
    Omega = diag(1, D)
    Sigma = D / N * Omega 
    Sigma_inv = solve(Sigma)
    alpha = rep(1, D) 
    mu_0 = rep(0, D)
    
    LIL_D[d_i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5)
    
    J = 5000
    N_approx = 10
    u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame 
    
    u_df_full = preprocess(u_samps, D)
    
    # skew_rpart = rpart(psi_u ~ ., u_df_full)
    # plot(skew_rpart)
    
    approx_skew = approx_lil(N_approx, D, u_df_full, J/N_approx)
    LIL_D_approx[d_i] = mean(approx_skew)
    
}



rbind(LIL_D, LIL_D_approx)


# ------------------------------------------------------------------------------












