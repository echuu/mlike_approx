

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)


setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")         




psi = function(u, prior) {
    
    # log prior
    logprior = sum(dunif(u, -5, 5, log = TRUE))
    
    # log-likelihood
    loglik =  -((u[prior$D]^2 + u[1]^2 - prior$a)^2 / prior$b)^2 - 
        sum(((u[1:(prior$D-1)]^2 + u[2:prior$D]^2 - prior$a) / prior$b)^2);
    
    return(- logprior - loglik)
}



B = 100
J         = B * 2000     # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 3000         # number of burn in draws
n_chains  = 8            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
a = 2
b = 1
D = 64
ring_dat = list(a = a, b = b, D = D) # for STAN sampler
params = list(a = a, b = b, D = D) # for STAN sampler

# (1) generate posterior samples -- should give us (J * N_approx) draws
ring_sampler = "C:/Users/ericc/mlike_approx/paper_simulations/results_ring/ring_stan.stan"
ring_fit = stan(file    =  ring_sampler, 
                data    =  ring_dat,
                iter    =  J_iter,
                warmup  =  burn_in,
                chains  =  8,                  
                control =  list(adapt_delta = 0.9999),  
                refresh = 0)     

ring_bridge = bridgesampling::bridge_sampler(ring_fit, silent = T)
ring_bridge$logml
LIL = -114.492

u_samp = rstan::extract(ring_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim

u_df_all = preprocess(u_post, D, params)

B = 100
xsplit     = rep(1:B, times = rep(J/B, B))
u_df_list = split.data.frame(u_df_all, xsplit)



u_df = u_df_list[[5]] # (J x D)
param = params


# move to FIXED_ALGO.R file 

J = J / B

# ------------------------------------------------------------------------------
source("setup.R")
set.seed(1)
hybrid = logml(D, u_df, J, param)
hybrid$param_out
if (any(is.na(hybrid))) {print(paste("error in iteration", b_i)); next;}
hybrid$all_approx
hybrid$wt_approx1 
hybrid$wt_approx2 
hybrid$wt_approx3
# ------------------------------------------------------------------------------


n_stage = 2
n_samps = 10
set.seed(1)
logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
prev_part = logml_approx$param_out$optimal_part # first stage partition
# psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)











