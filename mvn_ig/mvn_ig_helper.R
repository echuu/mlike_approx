

## log(determinant) function
# log_det = function(xmat) {
#     return(c(determinant(xmat, logarithm = T)$modulus))
# }


# lil() : compute the log marginal likelihood ----------------------------------
# y     : response vector (N x 1)
# X     : design matrix (N x d)
# d     : dimension of beta
# prior : list containing the prior parameter values
# post  : list containing the posterior paramter values
lil = function(y, X, prior, post, N = length(y), d = ncol(X)) {
    
    # u = (beta1,...,betad, sigmasq) \in R^(d+1)
    
    V_star = post$V_star
    mu_star = post$mu_star
    a_n = post$a_n
    b_n = post$b_n 
    
    V_beta = prior$V_beta
    mu_beta = prior$mu_beta
    a_0 = prior$a_0
    b_0 = prior$a_0
    
    log_py = a_0 * log(b_0) + lgamma(a_n) + 0.5 * log_det(V_star) - 
        0.5 * log_det(V_beta) - N / 2 * log(2 * pi) - lgamma(a_0) - 
        a_n * log(b_n)
    
    return(log_py)
    
} # end of lil() function


#### specify prior distribution density

## updated 1/12 -- is this used?
## log of the multivariate normal - inverse gamma density -- 
## log NIG(beta, sigmasq | mu_beta, V_beta, a, b)
log_mvnig = function(u, post, d = length(u) - 1) {
    
    mu_beta    = post$mu_star
    V_beta     = post$V_star
    V_beta_inv = post$V_star_inv
    a          = post$a_n
    b          = post$b_n
    
    # unname + unlist will remove the matrix/df structure and turn into vector,
    # -> easier handle, prevents input errors
    
    # print(d)
    
    beta = unname(unlist(u[1:d]))       # extract beta from the post sample
    sigmasq = unname(unlist(u[d + 1]))  # extract sigmasq from the post sample
    
    out = a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                           (beta - mu_beta))
    
    return(out %>% unlist() %>% unname())
}



## psi_true_mvn:     the true negative log posterior, not available in practice, but 
##                   we can evaluate it
## psi_mvn:          the negative log posterior as described in the notes, 
##                   = -loglik - logprior 

# done -- updated 1/12 -- don't think I ever use this? 
psi_true_mvn = function(u, post) {

    # p = length(u) - 1

    # beta = u[1:p]       # extract beta from the posterior sample
    # sigmasq = u[p + 1]  # extract sigmasq from the posterior sample

    # logpost = log_mvnig(u, post$mu_star, post$V_star, post$a_n, post$b_n)
    logpost = log_mvnig(u, post)

    return(-logpost %>% c()) # negative log posterior

    # return(p)
}

# checked + fixed -- 1/13
psi = function(u, prior) {
    
    y = prior$y
    X = prior$X
    
    d = length(u) - 1
    N = length(y)
    
    I_N = diag(1, N)       # (N x N) identity matrix
    
    # print(N)
    
    # extract prior parameters
    mu_beta    = prior$mu_beta
    V_beta     = prior$V_beta
    V_beta_inv = prior$V_beta_inv
    a          = prior$a_0
    b          = prior$b_0
    
    # extract beta, sigmasq from posterior sample -> used to evaluate the 
    # likelihood and the prior
    beta    = unname(unlist(u[1:d]))     # (p x 1) -- 1/13 : additional unname()
    sigmasq = unname(unlist(u[d+1]))     # (1 x 1) -- 1/13 : 
    
    
    # y_mu  = X %*% beta
    # y_var = sigmasq * I_N 
    
    # loglik = dmvnorm(c(y), mean = y_mu, sigma = y_var, log = T)
    
    # matches with:
    loglik = sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigmasq), log = T))
    
    # print(loglik)
    
    logprior = a * log(b) - d / 2 * log(2 * pi) - 
        0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                           (beta - mu_beta)) 
    
    # convert to scalar data type 
    logprior = logprior %>% c()
    
    # print(logprior)
    
    return(- loglik - logprior)
    
} # end of psi_mvn() function


# closed form expression for the gradient of psi_mvn()
# returns a D-dim vector containing the partial derivatives of each param
lambda = function(u, prior) {
    
    # extract priors
    y       = prior$y
    X       = prior$X
    mu_beta = prior$mu_beta
    V_beta  = prior$V_beta
    a       = prior$a_0
    b       = prior$b_0
    
    
    # closed form expression of the gradient
    
    N = length(y)
    p = length(u) - 1   # dimension of beta
    
    l_beta = numeric(p)      # derivative wrt beta
    l_sigmasq = numeric(1)   # derivative wrt sigmasq
    
    # evaluate gradient at these points
    beta = u[1:p]
    sigmasq = u[p + 1]
    
    l_beta = V_beta_inv %*% (beta - mu_beta) - t(X) %*% (y - X %*% beta)
    
    l_sigmasq = (N / 2 + a + p / 2 + 1) - 
        1 / sigmasq * (0.5 * t(y - X %*% beta) %*% (y - X %*% beta) + 
                           b + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                           (beta - mu_beta))
        
    return(1 / sigmasq * c(l_beta, l_sigmasq))
}


# preprocess = function(stan_fit, D, post, prior) {
#     
#     u_samp = rstan::extract(stan_fit, pars = c("beta", "sigmasq"), 
#                             permuted = TRUE)
#     
#     u_beta = u_samp$beta %>% data.frame()
#     u_sigmasq = u_samp$sigmasq
#     
#     u_post = data.frame(beta_post = u_beta, sigmasq_post = u_sigmasq)
#     
#     #### 2/10 update --- switch over to the psi that we should actually be using 
#     # psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)
#     psi_u = apply(u_post, 1, psi, prior = prior) %>% unname() # (J x 1)
#     
#     # (1.2) construct u_df -- this will require some automation for colnames
#     u_df_names = character(D + 1)
#     for (d in 1:D) {
#         u_df_names[d] = paste("u", d, sep = '')
#     }
#     u_df_names[D + 1] = "psi_u"
#     
#     # populate u_df
#     u_df = cbind(u_post, psi_u) # (J * N_approx) x (D + 1)
#     
#     # rename columns (needed since these are referenced explicitly in partition.R)
#     names(u_df) = u_df_names
#     
#     
#     return(u_df)
# 
# }


preprocess_resample = function(u_post, D, prior) {
    
    psi_u = apply(u_post, 1, psi, prior = prior) %>% unname() # (J x 1)
    
    # (1.2) construct u_df -- this will require some automation for colnames
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(u_post, psi_u) # (J * N_approx) x (D + 1)
    
    # rename columns (needed since these are referenced explicitly in partition.R)
    names(u_df) = u_df_names
    
    
    return(u_df)
}









