library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library("ggpmisc")
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)


setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")         


# pd.DataFrame(rec.get().logp).to_csv("ring_psi.csv")
# pd.DataFrame(rec.get().samples).to_csv("ring.csv")


setwd("C:/Users/ericc/Dropbox/eric chuu research/aistats/rdata_files")
ring_mat = read.csv("ring.csv")[,-1] %>% as.matrix
ring = read.csv("ring.csv")[,-1] %>% data.frame()
ring_psi = read.csv("ring_psi.csv")[,-1]


a = 2
b = 1
D = 20
prior_const = sum(dunif(rep(0, D), -5, 5, log = TRUE))

psi = function(u, prior) {
    
    # log prior
    # logprior = sum(dunif(u, -5, 5, log = TRUE))
    
    # log-likelihood
    # loglik =  -((u[prior$D]^2 + u[1]^2 - prior$a)^2 / prior$b)^2 -
    #     sum(((u[1:(prior$D-1)]^2 + u[2:prior$D]^2 - prior$a)^2 / prior$b)^2);
    
    # return(- logprior - loglik)
    u2 = u^2
    u2p = c(u2[D], u2, u2[1])
    sum((u2p[1:D] + u2p[2:(D+1)] - a)^2 / b) - prior_const
    
}


ring_psi %>% head # need to get my psi to match this calculation
# u = ring[1,1:D] %>% unname %>% unlist
# 
# ring_dat = list(a = a, b = b, D = D) # for STAN sampler
params = list(a = a, b = b, D = D) # for STAN sampler
# 
# 
# u2 = u^2
# u2p = c(u2[D], u2, u2[1])
# 
# u2p %>% length
# u2p[1:D] %>% length
# u2p[2:(D+1)] %>% length
# 
# sum((u2p[1:D] + u2p[2:(D+1)] - a)^2 / b) - sum(dunif(u, -5, 5, log = TRUE))

u_df_all = preprocess(ring, D, params)

u_df_all$psi_u %>% head

J = nrow(u_df_all)
hybrid = hybrid_ml(D, u_df_all[sample(1:J, 1000),], 1000, params)
hybrid$zhat


# bridge example
library(bridgesampling)
log_density = function(u, data) {
    -psi(u, data)
}

lb <- rep(-5, D)
ub <- rep(5, D)
u_samp = as.matrix(ring[sample(1:J, 1000),1:D])
colnames(u_samp) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridge_sampler(samples = u_samp,
                               log_posterior = log_density,
                               data = params, lb = lb, ub = ub, silent = TRUE)
bridge_result$logml


B = 100
J = 1000
xsplit    = rep(1:B, times = rep(J, B))
u_df_list = split.data.frame(u_df_all, xsplit)




i = 1

# LIL = -28.618178
#  LIL = -42.93973
LIL = -35.78178

set.seed(123)
bridge = numeric(B)
hyb = numeric(B)
for (i in 1:B) {
    
    # u_df = u_df_all[sample(1:nrow(u_df_all), 1500),]
    u_df = u_df_list[[i]]
    
    lb <- rep(-5, D)
    ub <- rep(5, D)
    u_samp = as.matrix(u_df[,1:D])
    colnames(u_samp) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(u_samp)
    
    bridge_result = bridge_sampler(samples = u_samp,
                                   log_posterior = log_density, data = params, 
                                   lb = lb, ub = ub, silent = TRUE)
    bridge_result$logml
    bridge[i] = bridge_result$logml
    
    
    #### (1) hybrid estimator
    hybrid = hybrid_ml(D, u_df, J, prior)
    
    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    hyb[i] = hybrid$zhat
    
    # avg_hyb = mean(hyb[hyb!=0])
    print(paste("iter ", i, ': ',
                "hybrid = ", round(mean(hyb[hyb!=0]), 3),
                '; ', "mae = ", round(mean((LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ',
                "mae = ", round(mean((LIL - bridge[bridge!=0])), 4),
                # "came = ", round(mean(came[came!=0]), 3), '; ', 
                # "mae = ", round(mean(abs(LIL - came[came!=0])), 4),
                sep = '')) 
}


approx = data.frame(LIL, hyb = hyb[hyb!=0], bridge = bridge[bridge!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans((LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)







B = 100
J         = B * 5000     # number of MC samples per approximation
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


u_samp = rstan::extract(ring_fit, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
u_post %>% dim
u_df = preprocess(u_post, D, params)
# getwd()
# save.image()

hml_approx = hml_const(N_approx, D, u_df_list[[3]], J/B, params)
hml_approx$const_vec
hybrid = hybrid_ml(D, u_df_list[[7]], J/B, params)
hybrid$zhat

n_samps = 10
# og_part = hml_approx$param_out %>%
#     dplyr::select(-c(psi_choice, logQ_cstar))
# ss_part = fit_resid(og_part, D, n_samps, params)
# ts_part = fit_resid(ss_part, D, n_samps / 2, params)
# log_sum_exp(unlist(compute_expterms(ss_part, D)))
# log_sum_exp(unlist(compute_expterms(ts_part, D)))
# mean(c(hml_approx$const_vec, 
#        log_sum_exp(unlist(compute_expterms(ss_part, D))),
#        log_sum_exp(unlist(compute_expterms(ts_part, D)))))


# ------------------------------------------------------------------------------
B = 50
# J / 1000

xsplit     = rep(1:B, times = rep(J/B, B))
u_df_list = split.data.frame(u_df, xsplit)


# ### troubleshooting hml_const() function call ----------------------------------
# 
# u_rpart = rpart(psi_u ~ ., u_df[1:2000,])
# 
# 
# # (3.1) obtain the (data-defined) support for each of the parameters
# param_support = extractSupport(u_df, D) #
# 
# # (3.2) obtain the partition
# u_partition = extractPartition(u_rpart, param_support) 
# 
# # organize all data into single data frame --> ready for approximation
# param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R

### troubleshooting hml_const() function call ----------------------------------

i = 1
hyb     = numeric(B) # store hybrid estimator
set.seed(1)
for(i in 104:B) {
    
    # hml_approx = hml_const(N_approx, D, u_df_list[[i]], J/B, params)
    hml_approx = hybrid(u_df_list[[i]])
    if (any(is.na(hml_approx))) {print(paste("error in iteration", i)); next;}

    og_part = hml_approx$param_out %>%
        dplyr::select(-c(psi_choice, logQ_cstar))
    # ss_part = fit_resid(og_part, D, n_samps, params)
    # ts_part = fit_resid(ss_part, D, n_samps / 2, params)
    
    ss_part = resample_partition(og_part, n_samps)
    if (any(is.na(ss_part))) {print(paste("error in iteration", i)); next;}
    
    ts_part = resample_partition(ss_part, n_samps / 2)
    if (any(is.na(ts_part))) {print(paste("error in iteration", i)); next;}
    
    #### note: if any hit error, no values are filled ####
    
    hyb_fs[i] = hml_approx$const_vec
    hyb_ss[i] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    hyb_ts[i] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    hyb[i] = mean(c(hyb_fs[i], hyb_ss[i], hyb_ts[i]))

    print(paste("iter: ", i, 
                "hyb = ", round(mean(hyb[hyb != 0]), 4),
                "hyb3 = ", round(mean(hyb_ts[hyb_ts != 0]), 4)))
        
    # if (i %% 1 == 0) { 
    #     print(paste("iter: ", i, 
    #                 "hyb = ", round(mean(hyb[hyb != 0]), 4),
    #                 "hyb3 = ", round(mean(hyb_ts[hyb_ts != 0]), 4)))
    # }
    # i = i + 1
}

LIL = -114.492
approx = data.frame(hyb = hyb[hyb!=0], hyb_ss = hyb_ss[hyb_ss!= 0], 
                    hyb_ts = hyb_ts[hyb_ts!=0], LIL)

# approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)


set.seed(1)
bridge  = numeric(B) # bridge estimator
for (b_i in 1:B) {
    
    ring_fit = stan(file    =  ring_sampler, 
                    data    =  ring_dat,
                    iter    =  J_iter,
                    warmup  =  burn_in,
                    chains  =  8,                  
                    control =  list(adapt_delta = 0.9999),  
                    refresh = 0)     
    
    ring_bridge = bridgesampling::bridge_sampler(ring_fit, silent = T)
    ring_bridge$logml
    
    bridge[b_i] = ring_bridge$logml
    
    
    print(paste("iter: ", b_i, "bridge = ", round(mean(bridge[1:b_i]), 4))) 
    # if (b_i %% 1 == 0) { 
    #     print(paste("iter: ", b_i, 
    #                 "bridge = ", round(mean(bridge[1:b_i]), 4))) 
    # }
}


hybrid_res = hybrid(u_df[1:2000,])
is.na(hybrid_res)



resample_partition = function(fit_partition, n_resamp) {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully
            
            # message("This is the 'try' part")
            fit_resid(fit_partition, D, n_samps, params)
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("resample computation error"))
            # message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
            # NOTE:
            # Here goes everything that should be executed at the end,
            # regardless of success or error.
            # If you want more than one expression to be executed, then you 
            # need to wrap them in curly brackets ({...}); otherwise you could
            # just have written 'finally=<expression>' 
            # message(paste("Processed URL:", url))
            # message("Some other message at the end")
        }
    )    
    return(out)
}




hybrid <- function(u_df) {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully
            
            # message("This is the 'try' part")
            hml_const(N_approx, D, u_df, J/B, params)
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            # message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning=function(cond) {
            message("Here's the original warning message:")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
            # NOTE:
            # Here goes everything that should be executed at the end,
            # regardless of success or error.
            # If you want more than one expression to be executed, then you 
            # need to wrap them in curly brackets ({...}); otherwise you could
            # just have written 'finally=<expression>' 
            # message(paste("Processed URL:", url))
            # message("Some other message at the end")
        }
    )    
    return(out)
}










