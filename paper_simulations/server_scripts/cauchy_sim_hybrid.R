
setwd("/home/grad/ericchuu/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library(bridgesampling)
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)


D = 48
mu = 5
sigma = 1
cauchy_dat = list(mu = mu, sigma = sigma, D = D) # for STAN sampler
params = list(mu = mu, sigma = sigma, D = D)     # for hybrid algo
LIL = -254.627
n_samps = 10



# ------------------------------------------------------------------------------

resample_partition = function(fit_partition, n_resamp) {
    out <- tryCatch(
        {
            fit_resid(fit_partition, D, n_resamp, params)
            
        },
        error=function(cond) {
            message(paste("resample computation error"))
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            return(NULL)
        },
        finally={
            
        }
    )    
    return(out)
}
hybrid <- function(u_df) {
    out <- tryCatch(
        {
            
            hml_const(N_approx, D, u_df, J, params)
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            return(NULL)
        },
        finally={
            
        }
    )    
    return(out)
}

B         = 100          
J         = 2000         # number of MC samples per approximation
burn_in   = 5000         # number of burn in draws
n_chains  = 8            # number of markov chains to run
stan_seed = 123          # seed
J_iter = 1 / n_chains * B * J + burn_in 
cauchy_fit = stan(file    =  cauchy_sampler, 
                  data    =  cauchy_dat,
                  iter    =  J_iter,
                  warmup  =  burn_in,
                  chains  =  n_chains,                  
                  control =  list(adapt_delta = 0.999),  
                  refresh = 0)    
u_samp = rstan::extract(cauchy_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim
u_df = preprocess(u_post, D, params)

xsplit     = rep(1:B, times = rep(J, B))
u_df_list = split.data.frame(u_df, xsplit)

# check dimension of list elements
u_df_list[[1]] %>% dim # should be (J x D)

hyb     = numeric(B) # store hybrid estimator
hyb_fs  = numeric(B) # store fs estimator
hyb_ss  = numeric(B) # store ss estimator
hyb_ts  = numeric(B) # store ts estimator
for(i in 1:B) {
    
    # hml_approx = hml_const(N_approx, D, u_df_list[[i]], J/B, params)
    # og_part = hml_approx$param_out %>%
    #     dplyr::select(-c(psi_choice, logQ_cstar))
    # ss_part = fit_resid(og_part, D, n_samps, params)
    # ts_part = fit_resid(ss_part, D, n_samps / 2, params)
    
    ## try catch all the algo calls to avoid crashing
    hml_approx = hybrid(u_df_list[[i]])
    if (any(is.na(hml_approx))) {print(paste("error in iteration", i)); next;}
    
    og_part = hml_approx$param_out %>% dplyr::select(-c(psi_choice, logQ_cstar))
    ss_part = resample_partition(og_part, n_samps)
    if (any(is.na(ss_part))) {print(paste("error in iteration", i)); next;}
    
    ts_part = resample_partition(ss_part, n_samps / 2)
    if (any(is.na(ts_part))) {print(paste("error in iteration", i)); next;}
    
    hyb_fs[i] = hml_approx$const_vec
    hyb_ss[i] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    hyb_ts[i] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    hyb[i] = mean(c(hyb_fs[i], hyb_ss[i], hyb_ts[i]))
    
    print(paste("iter: ", i, 
                "hyb = ", round(mean(hyb[hyb!=0]), 4),
                "hyb3 = ", round(mean(hyb_ts[hyb_ts!= 0]), 4)))
}

approx = data.frame(hyb = hyb[hyb!=0], hyb_ss = hyb_ss[hyb_ss!= 0], 
                    hyb_ts = hyb_ts[hyb_ts!=0], LIL)

# approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)

save.image()
