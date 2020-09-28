

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")  
library(rstan)

rstan_options(auto_write = TRUE)


# sample data
a = 1
b = 0.5
d = 16
params = list(a = a, b = b, D = d)
psi = function(u, prior) {
    
    # log prior
    logprior = dunif(u[1], -4, 4, log = TRUE) + 
        sum(dunif(u[-1], -30, 30, log = TRUE))
    
    # log-likelihood
    loglik = dnorm(u[1], 0, prior$a, log = T) + 
        sum(dnorm(u[-1], 0, sqrt(exp(2 * prior$b * u[1])), log = T))
    
    return(- logprior - loglik)
}

B = 100
bridge  = numeric(B) # bridge estimator (normal)

J         = 5000 * B     # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 5000         # number of burn in draws
n_chains  = 8            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
D = 16
funnel_dat = list(a = a, b = b, D = D) # for STAN sampler
params = list(a = a, b = b, D = D)

# (1) generate posterior samples -- should give us (J * N_approx) draws
funnel_sampler = "C:/Users/ericc/mlike_approx/paper_simulations/results_funnel/funnel_stan.stan"
funnel_fit = stan(file    =  funnel_sampler, 
                  data    =  funnel_dat,
                  iter    =  J_iter,
                  warmup  =  burn_in,
                  chains  =  n_chains,                  
                  control =  list(adapt_delta = 0.999),  
                  refresh = 0)     
u_samp = rstan::extract(funnel_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim
u_df = preprocess(u_post, D, params)

xsplit     = rep(1:B, times = rep(J/B, B))
u_df_list = split.data.frame(u_df, xsplit)

LIL = -63.4988

B = 10 # number of replications
hyb = numeric(B)
set.seed(1)
for (i in 1:B) {
    
    hybrid = hybrid_ml(D, u_df_list[[i]], J/B, params)
    
    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    
    hyb[i] = hybrid$zhat
    
    avg_hyb = mean(hyb[hyb!=0])
    print(paste("iter ", i, ': ',
                "hybrid = ", round(avg_hyb, 3), '; ',
                "ae = ", LIL - avg_hyb,
                sep = '')) 
}

approx = data.frame(LIL, hyb = hyb[hyb!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)

saveRDS(list(J = J, D = D, approx_df = approx), 
        file = 'funnel_d5.RData')
funnel_d16 = readRDS('funnel_d5.RData')






