
setwd("C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files")
funnel_mat = read.csv("funnel.csv")[,-1] %>% as.matrix
funnel = read.csv("funnel.csv")[,-1] %>% data.frame()
funnel_psi = read.csv("funnel_loglike.csv")[,-1]

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")

# sample data
a = 1
b = 0.5
D = 10
params = list(a = a, b = b, D = D)
psi = function(u, prior) {
    
    # log prior
    logprior = dunif(u[1], -4, 4, log = TRUE) + 
        sum(dunif(u[-1], -30, 30, log = TRUE))
    
    # log-likelihood
    loglik = dnorm(u[1], 0, prior$a, log = T) + 
        sum(dnorm(u[-1], 0, sqrt(exp(2 * prior$b * u[1])), log = T))
    
    return(- logprior - loglik)
}


u_df = preprocess(funnel, D, params)
J = u_df %>% nrow

u_df$psi_u %>% head
funnel_psi %>% head

hybrid = hybrid_ml(D, u_df[sample(1:J, 100),], 100, params)
hybrid$zhat


log_density = function(u, data) {
    -psi(u, data)
}

lb <- c(-4, rep(-30, D-1))
ub <- c(4, rep(30, D-1))
colnames(funnel_mat) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(funnel_mat)

library(bridgesampling)
bridge_result = bridge_sampler(samples = funnel_mat,
                               log_posterior = log_density,
                               data = params, lb = lb, ub = ub, 
                               silent = TRUE)
bridge_result$logml


log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])




# STAN SETTINGS ----------------------------------------------------------------
B = 10
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

funnel_bridge = bridgesampling::bridge_sampler(funnel_fit, silent = FALSE)
funnel_bridge$logml
funnel_bridge$niter

u_samp = rstan::extract(funnel_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim
u_df = preprocess(u_post, D, params)
hml_approx = hml_const(N_approx, D, u_df, J, params)
hml_approx$const_vec
funnel_bridge$logml

n_samps = 10
og_part = hml_approx$param_out %>%
    dplyr::select(-c(psi_choice, logQ_cstar))
ss_part = fit_resid(og_part, D, n_samps, params)
ts_part = fit_resid(ss_part, D, n_samps / 2, params)

log_sum_exp(unlist(compute_expterms(ss_part, D)))
log_sum_exp(unlist(compute_expterms(ts_part, D)))

mean(c(hml_approx$const_vec, 
       log_sum_exp(unlist(compute_expterms(ss_part, D))),
       log_sum_exp(unlist(compute_expterms(ts_part, D)))))




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

i = 1
hyb     = numeric(B) # store hybrid estimator
hyb_fs  = numeric(B) # store fs estimator
hyb_ss  = numeric(B) # store ss estimator
hyb_ts  = numeric(B) # store ts estimator
for(i in 14:B) {
    
    hml_approx = hml_const(N_approx, D, u_df_list[[i]], J/B, params)
    og_part = hml_approx$param_out %>%
        dplyr::select(-c(psi_choice, logQ_cstar))
    ss_part = fit_resid(og_part, D, n_samps, params)
    ts_part = fit_resid(ss_part, D, n_samps / 2, params)
    
    hyb_fs[i] = hml_approx$const_vec
    hyb_ss[i] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    hyb_ts[i] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    hyb = mean(c(hyb_fs[i], hyb_ss[i], hyb_ts[i]))
    
    if (i %% 10 == 0) { 
        print(paste("iter: ", i, 
                    "hyb = ", round(mean(hyb[1:i]), 4)))
    }
    i = i + 1
}



set.seed(1)
for (b_i in 1:B) {
    
    funnel_fit = stan(file    =  funnel_sampler, 
                      data    =  funnel_dat,
                      iter    =  J_iter,
                      warmup  =  burn_in,
                      chains  =  n_chains,                  
                      control =  list(adapt_delta = 0.999),  
                      refresh = 0)     
    
    # y = rstan::extract(funnel_fit, pars = c("u"), permuted = TRUE)
    funnel_bridge = bridgesampling::bridge_sampler(funnel_fit, silent = TRUE)
    bridge[b_i] = funnel_bridge$logml

    if (b_i %% 10 == 0) { 
        print(paste("iter: ", b_i, 
                    # "hyb = ", round(mean(bridge[1:b_i]), 3),
                    # "came = ", round(mean(came[1:b_i]), 3),
                    "bridge = ", round(mean(bridge[1:b_i]), 3))) 
    }
}


approx = data.frame(bridge, LIL)

# approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)





