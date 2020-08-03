


log_density = function(u, data) {
    -psi(u, data)
}




## (2) obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post

# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)
u_samp = as.matrix(post_samps)
colnames(u_samp) = names(u_df)[1:D_u]
# prepare bridge_sampler input()
lb = rep(-Inf, D_u)
ub = rep(Inf, D_u)

diag_ind = getDiagIndex(D, D_u) # obtain column index of the diagonal entries
lb[diag_ind] = 0                # diagonal entries are positive
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
                                               log_posterior = log_density,
                                               data = param_list, 
                                               lb = lb, ub = ub, 
                                               silent = TRUE)
bridge_result$logml
true_logml

B = 100
bridge = numeric(B)
set.seed(1)
for (b_i in 1:B) {
    
    postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post
    
    # these are the posterior samples stored as vectors (the lower cholesky factors
    # have been collapsed into D_u dimensional vectors)
    post_samps = postIW$post_samps                 # (J x D_u)
    
    
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
    # # set.seed(1)
    # ss_part = fit_resid(og_part, D_u, n_samps, param_list)
    # ts_part = fit_resid(ss_part, D_u, n_samps / 2, param_list)
    
    # hyb_fs[b_i] = hml_approx$const_vec
    # hyb_ss[b_i] = log_sum_exp(unlist(compute_expterms(ss_part, D_u)))
    # hyb_ts[b_i] = log_sum_exp(unlist(compute_expterms(ts_part, D_u)))
    
    u_samp = as.matrix(post_samps)
    colnames(u_samp) = names(u_df)[1:D_u]
    
    bridge_result = bridgesampling::bridge_sampler(samples = u_samp, 
                                                   log_posterior = log_density,
                                                   data = param_list, 
                                                   lb = lb, ub = ub, 
                                                   silent = TRUE)
    bridge[b_i] = bridge_result$logml
    
    if (b_i %% 10 == 0) { 
        print(paste("iter: ", b_i, ' - ',
                    round(mean(bridge[1:b_i]), 3), 
                    sep = '')) 
    }
    
}

LIL = true_logml
approx = data.frame(bridge)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)









