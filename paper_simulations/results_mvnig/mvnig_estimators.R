

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

A_mu = c(min(mu_post), max(mu_post))
A_sigmasq = c(min(sigma_sq_post), max(sigma_sq_post))

# draw from the importance density N-IG
sigmasq_mean = mean(sigma_sq_post)
sigmasq_var = var(sigma_sq_post)

# compute shape/scale parameters from posterior mean/var of sigmasq
r_imp = sigmasq_mean^2/sigmasq_var + 2
s_imp = sigmasq_mean * (1 + sigmasq_mean^2 / sigmasq_var)

mu_s      = rnorm(J, mean(mu_post), sqrt(mean(sigma_sq_post) / w_n))
sigmasq_s = MCMCpack::rinvgamma(J, shape = r_imp, scale = s_imp)

# compute 1/s(theta) -- (K x 1) vector of evaluated densities
s_theta = dnorm(mu_s, m_n, sqrt(sigma_sq / w_n)) * 
    MCMCpack::dinvgamma(sigmasq_s, shape = r_n / 2, scale = s_n / 2)

# compute prior density
p_theta = dnorm(mu_s, m_0, sqrt(sigma_sq / w_0)) * 
    MCMCpack::dinvgamma(sigmasq_s, shape = r_0 / 2, scale = s_0 / 2)




#### CORRECTED ARITHMETIC MEAN ESTIMATOR ---------------------------------------
came_approx = function(u_df, hml_approx, prior, post, J, D) {
    
    A_supp = hml_approx$param_support # data defined support (D x 2)
    
    ## draw from importance distribution
    post_mean = unlist(unname(colMeans(u_df[,1:D])))
    
    beta_mean    = post_mean[1:p]
    sigmasq_mean = post_mean[D]
    sigmasq_var  = var(u_df[,D])
    r_imp = sigmasq_mean^2/sigmasq_var + 2
    s_imp = sigmasq_mean * (1 + sigmasq_mean^2 / sigmasq_var)
    
    
    beta_imp = rmvnorm(J, mean = beta_mean, 
                       sigma = sigmasq_mean * post$V_star) %>% data.frame 
    sigmasq_imp = MCMCpack::rinvgamma(J, shape = r_imp, scale = s_imp)
    
    imp_samp = data.frame(beta_imp, sigmasq_imp)
    u_df_imp = preprocess(imp_samp, D, prior)
    
    # to be updated
    log_s_theta = unname(dmvnorm(beta_imp, mean = beta_mean, 
                                 sigma = sigmasq_mean * post$V_star, 
                                 log = TRUE)) + 
        log(MCMCpack::dinvgamma(sigmasq_imp, shape = r_imp, scale = s_imp))
    
    include_d = rep(TRUE, J)
    for (d in 1:D) {
        include_d = include_d & 
            (imp_samp[,d] >= A_supp[d,1]) & (imp_samp[,d] <= A_supp[d,2])
    }
    
    Jp = sum(include_d)
    
    came_j = -log(J) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
    came_jp = -log(Jp) + log_sum_exp((-u_df_imp$psi_u - log_s_theta)[include_d])
    
    return(c(came_j, came_jp))
}



#### HARMONIC MEAN ESTIMATOR ---------------------------------------------------

# reg_lik() function -----------------------------------------------------------
# reformulate the likelihood so that it is of the form exp(x) so that 
# we can take advantage of the log-sum-exp trick (below); this function
# returns the part of the reformulated likelihood that is in the exponential
reg_lik = function(u_df, prior, J, D, N) {
    
    # lik_inv = numeric(J) # store likelihood computed for MCMC samples
    tmp_j   = numeric(J) # store quantity that is passed into log(sum(exp(x)))
    p = D - 1
    
    sigmasq_samp = sigmasq
    
    for (j in 1:J) {
        
        beta_samp    = unname(unlist(u_df[j, 1:p]))
        # uncomment line below for NIG case
        sigmasq_samp = unname(unlist(u_df[j, D]))
        
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





