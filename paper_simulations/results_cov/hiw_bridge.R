

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/covariance")
source("HIW_helper.R")

log_density = function(u, data) {
    -psi(u, data)
}



post_samps = postIW$post_samps                 # (J x D_u)

post_samps %>% head

u_samp = as.matrix(post_samps)
colnames(u_samp) = names(u_df)[1:D_u]
# prepare bridge_sampler input()
lb = rep(-Inf, D_u)
ub = rep(Inf, D_u)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
                                               log_posterior = log_density,
                                               data = params, 
                                               lb = lb, ub = ub, 
                                               silent = TRUE)
bridge_result$logml


logmarginal(Y, testG, b, V, S)
gnorm_approx
hml_approx$const_vec


B = 100
J = 1000
hyb     = numeric(B) # store harmonic mean estiator
hyb_fs  = numeric(B) # store harmonic mean estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator

bridge  = numeric(B)
set.seed(1)
for (b_i in 1:B) {
    
    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)  
    post_samps = postIW$post_samps                 # (J x D_u)
    
    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
    hml_approx = hml_const(1, D_u, u_df, J, params)

    hyb_fs[b_i] = hml_approx$const_vec
    
    # u_df stores the posterior samples row-wise so that the first D_u columns 
    # store the lower cholesky factors in vector form, and the last column is
    # the function evaluate psi(u), so u \in R^(D_u), and psi(u) \in R
    # u_df = preprocess(post_samps, D_u, param_list) # J x (D_u + 1)
    # 
    # ## (3b) compute approximation
    # hml_approx = hml_const(1, D_u, u_df, J, param_list)
    # og_part = hml_approx$param_out %>%
    #     dplyr::select(-c(psi_choice, logQ_cstar))
    # 
    # # # set.seed(1)
    # ss_part = fit_resid(og_part, D_u, n_samps, params)
    # ts_part = fit_resid(ss_part, D_u, n_samps / 2, params)
    
    # hyb_fs[b_i] = hml_approx$const_vec
    # hyb_ss[b_i] = log_sum_exp(unlist(compute_expterms(ss_part, D_u)))
    # hyb_ts[b_i] = log_sum_exp(unlist(compute_expterms(ts_part, D_u)))
    
    # u_samp = as.matrix(post_samps)
    # colnames(u_samp) = names(u_df)[1:D_u]
    # 
    # bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
    #                                                log_posterior = log_density,
    #                                                data = params,
    #                                                lb = lb, ub = ub,
    #                                                silent = TRUE)
    # bridge[b_i] = bridge_result$logml
    
    if (b_i %% 1 == 0) { 
        print(paste("iter: ", b_i, ' - ',
                    # round(mean(bridge[1:b_i]), 3), ', ',
                    round(mean(hyb_fs[1:b_i]), 4),
                    sep = '')) 
    }
    
}



LIL = logmarginal(Y, testG, b, V, S)
approx = data.frame(hyb)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)



LIL = logmarginal(Y, testG, b, V, S)
approx = data.frame(bridge, hyb = hyb_fs)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)




