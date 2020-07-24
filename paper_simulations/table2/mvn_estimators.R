


#### CORRECTED ARITHMETIC MEAN ESTIMATOR ---------------------------------------
came_approx = function(u_df, hml_approx, prior, post, J, D) {
    
    A_beta = hml_approx$param_support # data defined support (D x 2)
    post_mean = unlist(unname(colMeans(u_df[,1:D])))
    
    imp_samp = rmvnorm(J, mean = post_mean, sigma = post$Q_beta_inv) %>% 
        data.frame 
    
    u_df_imp = preprocess(imp_samp, D, prior)
    
    log_s_theta = unname(dmvnorm(imp_samp, mean = post_mean, 
                                 sigma = post$Q_beta_inv, 
                                 log = TRUE))
    
    include_d = rep(TRUE, D)
    for (d in 1:D) {
        include_d = include_d & 
            (imp_samp[,d] >= A_beta[d,1]) & (imp_samp[,d] <= A_beta[d,2])
    }
    
    -log(J) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
}



#### HARMONIC MEAN ESTIMATOR ---------------------------------------------------

# reg_lik() function -----------------------------------------------------------
# reformulate the likelihood so that it is of the form exp(x) so that 
# we can take advantage of the log-sum-exp trick (below); this function
# returns the part of the reformulated likelihood that is in the exponential
reg_lik = function(u_df, prior, J, D, N) {
    
    # lik_inv = numeric(J) # store likelihood computed for MCMC samples
    tmp_j   = numeric(J) # store quantity that is passed into log(sum(exp(x)))
    p = D
    
    sigmasq_samp = sigmasq
    
    for (j in 1:J) {
        
        beta_samp    = unname(unlist(u_df[j, 1:p]))
        # uncomment line below for NIG case
        # sigmasq_samp = unname(unlist(u_df[j, p+1]))
        
        # lik_inv[j] = 1 / prod(dnorm(data$y, data$X %*% beta_samp, 
        #                         sqrt(sigmasq_samp)))
        
        tmp_j[j] = 1 / (2 * sigmasq_samp) * 
            sum((prior$y - prior$X %*% beta_samp)^2) + N / 2 * log(sigmasq_samp)
        
    } # end of loop iterating over MCMC samples
    
    
    return(tmp_j)
    
} # end of reg_lik() function --------------------------------------------------



# hme() function ---------------------------------------------------------------
# harmonic mean estimator -- this is written specifically for the MVN-IG
# example, since the likelihood is of a form such that we can take advantage
# of the log-sum-exp trick to stabilize the calculation of estimator
hme_approx = function(u_df, prior, J, D, N) {
    
    # harmonic mean estimator requires calculating the likelihood given
    # each of the J parameters (that are sampled via MCMC)
    
    # in order to generalize this function, each model that wants to take
    # advantage of the hme estimator should provide
    # (1) likelihood function
    # (2) parameter extraction that only requires u_df input
    
    # lik_inv = reg_lik(u_df, data, J, D)
    # hme_estimate = log(J) - log(sum(lik_inv))
    
    # log_sum_exp version of 
    tmp_j = reg_lik(u_df, prior, J, D, N)
    hme_estimate = log(J) - N / 2 * log(2 * pi) - log_sum_exp(tmp_j)
    
    return(hme_estimate)
    
} # end of hme() function ------------------------------------------------------
