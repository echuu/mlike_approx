

library(mvtnorm)

set.seed(1)

D = 3
N = 50

b_0 = rep(0, D)        # mu_beta
V_0 = diag(1, D)       # V_beta
r_0 = 2 / 2            # a
s_0 = 1 / 2            # b 

beta = c(5, 1, -2)

I_N = diag(1, N)
I_D = diag(1, D)

X = matrix(rnorm(N * D), N, D)

eps = rmvnorm(N, mean = 0)

y = X %*% beta + eps # (N x 1) response vector


# TODO: compute true LIL


#### specify prior distribution density

## log(determinant) function
log_det = function(xmat) {
    return(c(determinant(xmat, logarithm = T)$modulus))
}


## log of the multivariate normal - inverse gamma density -- 
## log NIG(beta, sigmasq | mu_beta, V_beta, a, b)
log_mvnig = function(u, mu_beta, V_beta, a, b, d = length(u) - 1) {
    
    # p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
    #     sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
    #                                                solve(V_beta) %*% 
    #                                                (beta - mu_beta)))
    
    beta = u[1:d]
    sigmasq = u[d+1]
    
    a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% solve(V_beta) %*% 
                           (beta - mu_beta))
}

psi_mvn = function(u, y, mu_beta, V_beta, a, b, 
                   n = length(y), d = length(u) - 1) {

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



## psi_true:     the true negative log posterior, not available in practice, but 
##               we can evaluate it
## psi:          the negative log posterior as described in the notes, 
##               = -loglik - logprior 
## psi_tilde:    approximation of psi as described in the notes
##               = c_k + lambda_k'u

psi_true_mvn = function(u, mu_star, V_star, a_n, b_n) {
    
    # -log (true) posterior
    # -log(MVN-inverse gamma density)
    
    beta = u[1:d]
    sigmasq = u[d+1]
    
    logpost = log_mvnig(u, mu_star, V_star, a_n, b_n, d = length(u) - 1)
    
    return(-logpost) # negative log posterior
    
} # end of psi_true_mvn() function



# lambda = grad(psi)
# u is the "representative point" of each partition from the tree output
# mu_beta, V_beta, a, b are the PRIOR PARAMTERS since psi_mvn evaluates the
# log-likelihood and the log-prior
lambda_mvn = function(u, y, mu_beta, V_beta, a, b) {
    
    
    # don't want to compute this by hand, so we use numerical differentiation
    # see compute_mlik_v1.R for example of grad() function
    # check: 1-d case matches the closed form that we calculated in nig_2d.R
    
    grad(psi_mvn, u, y = y, mu_beta = mu_beta, V_beta = V_beta, a = a, b = b)
    
} # end of lambda_mvn() function



# ------------------------------------------------------------------------------


library("numDeriv")
library('MCMCpack')   # for rinvgamma() function








# testing ----------------------------------------------------------------------



beta = c(5, 1, -2)
sigmasq = 4

u = c(beta, sigmasq)

log_mvnig(u, b_0, V_0, r_0, s_0)

grad(psi_mvn, u, y = y, mu_beta = b_0, V_beta = V_0, a = r_0, b = s_0)

psi_mvn(u, y = y, mu_beta = b_0, V_beta = V_0, a = r_0, b = s_0)



## TODO: scale the 2-d code so that we  get same results using refactored code





## TODO: check tree output for paramters of the form [beta, sigmasq],
#        where we rely on the order of the first d betas for later parts
#        of the algorithm, e.g., u[1:d], u[d+1]

## thought 1: I think before getting sent into the tree, we already pre-label
## the parameters as u1,...,ud so that we have a handle on every parameter



## TODO: try fitting tree for d' = 3, (beta1, beta2, sigmasq)
## thought 1: previous design just saw me manually creating u_df, i.e,
#  u1 = , u2 = , ... , up = . --- need a way to  avoid doing this, because i 
#  also then compute the closed form integral manually too -> O(d)










# ------------------------------------------------------------------------------














# old stuff --------------------------------------------------------------------


## log of the multivariate normal - inverse gamma density -- 
## log NIG(beta, sigmasq | mu_beta, V_beta, a, b)
log_mvnig = function(beta, sigmasq, mu_beta, V_beta, a, b, d = length(beta)) {
    
    # p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
    #     sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
    #                                                solve(V_beta) %*% 
    #                                                (beta - mu_beta)))
    
    a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% solve(V_beta) %*% 
                           (beta - mu_beta))
}



mvnig = function(beta, sigmasq, mu_beta, V_beta, a, b, d = length(beta)) {
    
    
    p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
        sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
                                                   solve(V_beta) %*% 
                                                   (beta - mu_beta)))
    
    return(p)
}



psi_mvn = function(u, y, n = length(y), d = length(u)) {
    
    loglik = dmvnorm(c(y), mean = X %*% u[1:(d-1)], sigma = u[d] * I_N, log = T)
    logprior = dmvnorm(u[1:(d-1)], mean = b_0, sigma = u[d] * I_D, log = T) + 
        log(dinvgamma(u[d], shape = r_0, scale = s_0))
    
    # print(loglik)
    # print(logprior)
    
    ## TODO: write out log(invgamma) density in closed form to prevent underflow
    
    
    -loglik - logprior
}




