


# pd.DataFrame(rec.get().logp).to_csv("funnel_loglike.csv")
# pd.DataFrame(rec.get().samples).to_csv("funnel.csv")

setwd("C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files")
banana_mat = read.csv("banana.csv")[,-1] %>% as.matrix
banana = read.csv("banana.csv")[,-1] %>% data.frame()
banana_psi = read.csv("banana_loglike.csv")[,-1]

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")



D = 32
Q = 0.01
banana_dat = list(Q = Q, D = D) # for STAN sampler

psi = function(u, prior) {
    
    # log prior
    logprior = sum(dunif(u, -15, 15, log = TRUE)) 
    
    # log-likelihood
    odd = seq(1, prior$D, by = 2)
    even = seq(2, prior$D, by = 2)
    loglik = -sum((u[odd]^2 - u[even])^2 / Q + (u[odd] - 1)^2)
    
    return(- logprior - loglik)
}
params = list(Q = Q, D = D)


u_df = banana %>% mutate(psi_u = -banana_psi)
names(u_df)[1:D] = paste('u', 1:D, sep = '_')
banana_psi %>% head

hybrid = hybrid_ml(D, u_df[sample(1:J, 1000),], 1000, banana_dat)
hybrid$zhat


B = 5
# STAN SETTINGS ----------------------------------------------------------------
J         = 2000 * B     # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 5000         # number of burn in draws
n_chains  = 8            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 




# (1) generate posterior samples -- should give us (J * N_approx) draws
banana_sampler = "C:/Users/ericc/mlike_approx/paper_simulations/results_banana/banana_stan.stan"
banana_fit = stan(file    =  banana_sampler, 
                  data    =  banana_dat,
                  iter    =  J_iter,
                  warmup  =  burn_in,
                  chains  =  8,                  
                  control =  list(adapt_delta = 0.9999),  
                  refresh = 0)     

banana_bridge = bridgesampling::bridge_sampler(banana_fit, silent = FALSE)
banana_bridge$logml


# hybrid approx here -----------------------------------------------------------
psi = function(u, prior) {
    
    # log prior
    logprior = sum(dunif(u, -15, 15, log = TRUE)) 

    # log-likelihood
    odd = seq(1, prior$D, by = 2)
    even = seq(2, prior$D, by = 2)
    loglik = -sum((u[odd]^2 - u[even])^2 / Q + (u[odd] - 1)^2)
    
    return(- logprior - loglik)
}
params = list(Q = Q, D = D)


u_samp = rstan::extract(banana_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim
u_df = preprocess(u_post, D, params)
hml_approx = hml_const(N_approx, D, u_df, J, params)
hml_approx$const_vec














set.seed(1)
bridge = numeric(B)
for (b_i in 1:B) {
    
    funnel_fit = stan(file    =  funnel_sampler, 
                      data    =  funnel_dat,
                      iter    =  J_iter,
                      warmup  =  burn_in,
                      chains  =  n_chains,                  
                      control =  list(adapt_delta = 0.999),  
                      refresh = 0)     
    funnel_bridge = bridgesampling::bridge_sampler(funnel_fit, silent = TRUE)
    bridge[b_i] = funnel_bridge$logml
    
    if (b_i %% 10 == 0) { 
        print(paste("iter: ", b_i, 
                    # "hyb = ", round(mean(bridge[1:b_i]), 3),
                    # "came = ", round(mean(came[1:b_i]), 3),
                    "bridge = ", round(mean(bridge[1:b_i]), 3))) 
    }
}

LIL = -63.4988
approx = data.frame(bridge, LIL)

# approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)














