

## log(determinant) function
log_det = function(xmat) {
    return(c(determinant(xmat, logarithm = T)$modulus))
}


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
        0.5 * log_det(V_star) - N / 2 * log(2 * pi) - lgamma(a_0) - 
        a_n * log(b_n)
    
    
    
    return(log_py)
    
} # end of lil() function


#### specify prior distribution density

## updated 1/12
## log of the multivariate normal - inverse gamma density -- 
## log NIG(beta, sigmasq | mu_beta, V_beta, a, b)
log_mvnig = function(u, post, d = length(u) - 1) {
    
    ## TODO: verify equality of values below
    ## TODO: put posterior params into post object
    
    # p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
    #     sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
    #                                                solve(V_beta) %*% 
    #                                                (beta - mu_beta)))
    
    mu_beta = post$mu_star
    V_beta  = post$V_star
    a       = post$a_n
    b       = post$b_n
    
    beta = unlist(u[1:d])  # extract beta from the posterior sample
    sigmasq = u[d + 1]     # extract sigmasq from the posterior sample
    
    # print(b)
    
    out = a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% solve(V_beta) %*% 
                           (beta - mu_beta))
    
    return(out %>% unlist() %>% unname())
}


## psi_true_mvn:     the true negative log posterior, not available in practice, but 
##                   we can evaluate it
## psi_mvn:          the negative log posterior as described in the notes, 
##                   = -loglik - logprior 

## psi_tilde:    approximation of psi as described in the notes
##               = c_k + lambda_k'u  (this is calculated in main loop)


# psi_true_mvn_0 = function(u, mu_star, V_star, a_n, b_n) {
    
    # -log (true) posterior
    # -log(MVN-inverse gamma density)
    
    ## TODO: verify this is calculating what's intended
    
#    beta = u[1:d]
#    sigmasq = u[d+1]
    
#    logpost = log_mvnig(u, mu_star, V_star, a_n, b_n, d = length(u) - 1)
    
    
    
#    return(-logpost) # negative log posterior
    
#} # end of psi_true_mvn() function


# updated 1/12
psi_true_mvn = function(u, post) {
    
    # p = length(u) - 1
    
    # beta = u[1:p]       # extract beta from the posterior sample
    # sigmasq = u[p + 1]  # extract sigmasq from the posterior sample
    
    # logpost = log_mvnig(u, post$mu_star, post$V_star, post$a_n, post$b_n)
    logpost = log_mvnig(u, post)
    
    return(-logpost) # negative log posterior
        
    # return(p)
}



# TODO: to be checked
psi_mvn = function(u, y, X, mu_beta, V_beta, a, b, 
                   n = length(y), d = length(u) - 1) {
    
    ## TODO: put *prior* params into prior object
    ## TODO: verify this is calculating what's intended
    
    beta = u[1:d]
    sigmasq = u[d+1]
    
    
    loglik = dmvnorm(c(y), mean = X %*% beta, sigma = sigmasq * I_N, log = T)
    
    logprior = a * log(b) - d / 2 * log(2 * pi) - 
        0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% solve(V_beta) %*% 
                           (beta - mu_beta))
    
    - loglik - logprior
} # end of psi_mvn() function


# TODO: to be checked
# lambda = grad(psi)
# u is the "representative point" of each partition from the tree output
# mu_beta, V_beta, a, b are the PRIOR PARAMTERS since psi_mvn evaluates the
# log-likelihood and the log-prior
lambda_mvn = function(u, y, mu_beta, V_beta, a, b) {
    
    
    # don't want to compute this by hand, so we use numerical differentiation
    # see compute_mlik_v1.R for example of grad() function
    
    ## TODO: closed form of the gradient is easy, check it by hand first
    
    grad(psi_mvn, u, y = y, mu_beta = mu_beta, V_beta = V_beta, a = a, b = b)
    
} # end of lambda_mvn() function











