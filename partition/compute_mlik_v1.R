
library('dplyr')      # data mangement
library('ggplot2')    # plotting
library('MCMCpack')   # for rinvgamma() function
library('rpart')      # for rpart() function
library('rpart.plot') # for rpart.rules() function
library('tidyr')


setwd("/home/eric/mlike_approx/partition")

psi_true = function(mu, sigma_sq, m_n, w_n, r_n, s_n) {
    
    # p (mu, sigma_sq | y ) 
    #    = N (mu | m_n, sigma_sq / w_n) * IG (sigma_sq | r_n / 2, s_n / 2)
    #    = NIG (mu, sigma_sq | )
    
    log_mu_pdf   = dnorm(mu, m_n, sqrt(sigma_sq / w_n), log = T)
    log_sigma_sq = log(MCMCpack::dinvgamma(sigma_sq, r_n / 2, s_n / 2))
    
    log_p_mu_sigmasq = log_mu_pdf + log_sigma_sq
    
    return(-log_p_mu_sigmasq)
}

# generate samples from the posterior probability to form the HME estimator
J = 5000 # number of random draws used per estimate

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

# (2) for each u drawn from the posterior, evaluate psi(u) = psi(mu, sigma_sq)
# the actual value of psi that we will be using for this example - previous
# definition of psi was not exactly as described in the notes
psi_0 = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n) # (J x 1)

# u_df = data.frame(mu = mu_post, sigsq = sigma_sq_post, psi_u = psi_0) # (J x 3)

# input for paramPartition() MUST have parameter names u1, u2, ... up
u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_0)    # (J x 3)

# fit decision tree
### use rpart to fit partition
nig_rpart = rpart(psi_u ~ ., u_df)

# plot the tree (matches general structure of the tree returned from tree())
plot(nig_rpart)
text(nig_rpart, cex = 0.7)

source("partition.R")

### obtain partition
nig_support = rbind(c(min(mu_post), max(mu_post)),
                    c(min(sigma_sq_post), max(sigma_sq_post)))
nig_partition = paramPartition(nig_rpart, nig_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(nig_rpart, u_df, nig_partition)

# define psi_hat function -- using likelihood, prior
psi_hat = function(mu_star, sigmasq_star, y, m_0, w_0, r_0, s_0) {
    
    loglik = sum(dnorm(y, mu_star, sqrt(sigmasq_star), log = TRUE))
    logprior = dnorm(mu_star, m_0, sqrt(sigmasq_star / w_0), log = TRUE) *
        + log(MCMCpack::dinvgamma(sigmasq_star, r_0 / 2, s_0 / 2))
    
    return(-loglik - logprior)
}

# define lambda function
lambda = function(mu_star, sigmasq_star, y, m_0, w_0, r_0, s_0) {
    
    n = length(y)
    
    lambda1 = -1 / sigmasq_star * sum((y - mu_star)) + 
        w_0 / sigmasq_star * (mu_star - m_0)
    
    lambda2 = n / (2 * sigmasq_star)  - 
        1 / (2 * sigmasq_star^2) * sum((y - mu_star)^2) +
        (r_0 / 2 + 3 / 2) / sigmasq_star - 
        1 / (2 * sigmasq_star^2) * (w_0 * (mu_star - m_0)^2 + s_0)
    
    return(c(lambda1, lambda2))
    
}

# function uses lambda_k() and psi_true() to compute closed form integrals over
# each of the partitions
n_partitions = nrow(nig_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)
for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    c_k[k] = exp(-psi_hat(param_out[k,]$u1_star, 
                          param_out[k,]$u2_star,
                          y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    l_k = lambda(param_out[k,]$u1_star, param_out[k,]$u2_star,
                 y, m_0, w_0, r_0, s_0)
    
    # 1st param calculation
    p1 = -1/l_k[1] * exp(-l_k[1] * (param_out[k,]$u1_ub - param_out[k,]$u1_lb))
    
    # 2nd param calculation
    p2 = -1/l_k[2] * exp(-l_k[2] * (param_out[k,]$u2_ub - param_out[k,]$u2_lb))
    
    zhat[k] = c_k[k] * p1 * p2
}

log(sum(zhat))





### ----------------------------------------------------------------------------

library("numDeriv")


phi = function(u, y) {
    n = length(y)
    n / 2 * log(2 * pi) + n / 2 * log(u[2]) + 
        1 / (2 * u[2]) * sum((y - u[1])^2)
}

log_prior = function(u) {
    dnorm(u[1], mean = m_0, sd = sqrt(u[2] / w_0), log = T) + 
        log(dinvgamma(u[2], shape = r_0 / 2, scale = s_0 / 2))
}

grad_log_prior = function(u) { 
    
    grad_mu = - w_0 / u[2] * (u[1] - m_0)
    grad_sigmasq = - (r_0 / 2 + 3 / 2) / u[2] +
        1 / (2 * u[2]^2) * (w_0 * (u[1] - m_0)^2 + s_0)
    return(c(grad_mu, grad_sigmasq))
}


grad_phi = function(mu, sigmasq, y) {
    n = length(y)
    
    phi_mu = -1/sigmasq * sum(y - mu)
    phi_sigmasq = n / (2 * sigmasq)  - 1 / (2 * sigmasq^2) * sum((y - mu)^2)
    
    return(c(phi_mu, phi_sigmasq))
}


u_k = c(param_out[1,]$u1_star, param_out[1,]$u2_star)

grad_phi(u_k[1], u_k[2], y = y) - grad_log_prior(u_k)

grad(phi, u_k, y = y) - grad(log_prior, u_k)






