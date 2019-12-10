
## compute_mlik.R --------------------------------------------------------------
## (1) compute the marginal likelihood using existing techniques (MCMC)
## (2) approximate marginal likelihood using new method

library('dplyr')      # data mangement
library('ggplot2')    # plotting
library('MCMCpack')   # for rinvgamma() function
library('rpart')      # for rpart() function
library('rpart.plot') # for rpart.rules() function
library('tidyr')

setwd("/home/eric/mlike_approx/partition")

## example: use example for which we know the value of the marginal likelihood
## 

## we use example from the lenk paper which has the harmonic mean estimator
## plotted over the true value of the log-integrated-likelihood 
## log-integrated-likelihood ?= log marginal likelihood


## (0) generate data -----------------------------------------------------------

set.seed(123)

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N = 200

# generate 50 samples from N(mu, sigma_sq)
y = rnorm(N, mu, sqrt(sigma_sq))

ybar = mean(y)

# compute posterior parameters
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2

p_y = pi^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
    s_0^(r_0 / 2) / s_n^(r_n / 2)

LIL = log(p_y) # -113.143 (paper says -117, but difference arises from RNG)

# ------------------------------------------------------------------------------

## (1) compute the LIL using hme (DONE) ----------------------------------------

J = 1000 # number of random draws used per estimate
B = 1000 # number of batch estimators


# compute the estimate for the i-th batch
lil_hat = numeric(B) # store the log integrated likelihood for each batch
for (b in 1:B) {
    
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
    # (2) compute the harmonic mean estimator
    lik_j = numeric(J)
    for (j in 1:J) {
        # (2.1) compute the likelihood under the j-th posterior parameters
        lik_j[j] = 1 /  prod(dnorm(y, mu_post[j], sqrt(sigma_sq_post[j])))
    }
    
    # form the b-th HM estimator for the log integrated likelihood
    lil_hat[b] = 1 / (1 / J * sum(lik_j))
    
}

# (B x 3) -- batch number, lil estimate, actual lil value
hme_df = data.frame(mcmc = 1:B, hme = log(lil_hat), lil = LIL)

# mean of the hme estimate over B batches -> -115.548

# ------------------------------------------------------------------------------


## (2) given a set of points in the parameter space, approximate the LIL
## using method described in notes

## (2.1) generate samples from the *exact* posterior (one of our assumptions
## is that we can do this) 

J = 1000 # number of MCMC draws from posterior

# (2.1.1) sample from q(u_1 | - ) :=  mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (2.1.2) sample from q(u_2 | - ) := sigma_sq | y
sigma_sq_post = rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

# library("invgamma")


## (2.2) define psi(u), and evaluate psi(u) for each of the MCMC samples
# y        : (N x 1)
# mu       : (J x 1) J samples of mu from the posterior
# sigma_sq : (J x 1) J samples of sigma^2 from the posterior
# m_0      : (1 x 1) prior mean
# w_0      : (1 x 1) scale the prior variance
# r_0      : (1 x 1) prior shape
# s_0      : (1 x 1) prior scale
psi = function(y, mu, sigma_sq, m_0, w_0, r_0, s_0) {
    
    # for  y_1:N, evaluate the likelihood for j-th set of parameter values
    loglik = numeric(length(mu))
    
    for (j in 1:length(mu)) {
        loglik[j] = sum(dnorm(y, mu[j], sqrt(sigma_sq[j]), log = TRUE))
    }
    
    # following two computations are unused for now
    log_p_mu = dnorm(mu, m_0, sqrt(sigma_sq / w_0), log = TRUE)
    log_p_sigma_sq = dinvgamma(sigma_sq, r_0 / 2, s_0 / 2)
    
    return(loglik) 
} # end psi() function

psi_u = psi(y, mu_post, sigma_sq_post, m_0, w_0, r_0, s_0)


## (2.3) send these through the tree via rpart()

## before doing all the parameter extraction, representative point extraction,
## look at the partitioniong done from rpart, tree to see if the partitions are
## doing what we expect

# (2.3.1) form u_df that has the response variable psi_u, with u1,...,up
u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u) # (J x 3)

# (2.3.2) fit rpart()
nig_rpart = rpart(psi_u ~ ., u_df)

# plot the tree (matches general structure of the tree returned from tree())
plot(nig_rpart)
text(nig_rpart, cex = 0.7)


## (2.4) obtain the representative points of each partition
source("partition.R")
nig_support = rbind(c(-Inf, Inf), c(0, Inf))
nig_test = paramPartition(nig_rpart, nig_support)  # partition.R


## running code out of test_partition.R -- TODO: turn into function later
u_rpart = nig_rpart



## (2.5) obtain the approximation of the LIL using method described in notes
## TODO: this still needs to be implemeneted 


# ------------------------------------------------------------------------------

# (1) re-plot using exact negative log posterior for psi function

# returns negative log posterior : -log(p(u|y))
# usually only -log(p(y|u)) - log(p(u))
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
psi_u = psi_p(y, mu_post, sigma_sq_post, m_0, w_0, r_0, s_0) # (J x 1)

# the actual value of psi that we will be using for this example - previous
# definition of psi was not exactly as described in the notes
psi_0 = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n) # (J x 1)

u_df = data.frame(mu = mu_post, sigsq = sigma_sq_post, psi_u = psi_0) # (J x 3)

# input for paramPartition() MUST have parameter names u1, u2, ... up
u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_0)    # (J x 3)

# (3) fit decision tree
u_tree = tree(psi_u ~ ., u_df)

# (4) plot the partition over the parameter space
par(mfrow = c(1,1))
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.9, col = "green",
     xlab = 'mu', ylab = 'sigma_sq', main = 'N = 50, J = 1000')
partition.tree(u_tree, add = TRUE, cex = 0.01, ordvars = c("u1", "u2"))


## TODO ------------------------------------------------------------------------

source("partition.R")

# (2) function to evaluate lambda_dk

# gradient of psi -- return 2 dim vector
lambda_k = function(u, param) {
    
    lambda1 = 1 / u$sigma_sq * (u$mu - param$m_n)
    lambda2 = (param$r_n / 2 + 3 / 2) * 1 / u$sigma_sq - 
        1 / u$sigma_sq^2 * (param$s_n + (u$mu - param$m_n)^2)
    
    return(c(lambda1, lambda2))
}

### use rpart to fit partition
nig_rpart = rpart(psi_u ~ ., u_df)

# plot the tree (matches general structure of the tree returned from tree())
plot(nig_rpart)
text(nig_rpart, cex = 0.7)

### obtain partition
nig_support = rbind(c(min(mu_post), max(mu_post)),
                    c(min(sigma_sq_post), max(sigma_sq_post)))
nig_partition = paramPartition(nig_rpart, nig_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(nig_rpart, u_df, nig_partition)

### form approximation

# function uses lambda_k() and psi_true() to compute closed form integrals over
# each of the partitions
n_partitions = nrow(nig_partition)
c_k = numeric(n_partitions)
l_k = numeric(n_partitions)
zhat = numeric(n_partitions)
for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    c_k[k] = exp(-psi_true(param_out[k,]$u1_star, 
                           param_out[k,]$u2_star,
                           m_n, w_n, r_n, s_n)) # (1 x 1)
    
    lambda = lambda_k(u = list(mu = param_out[k,]$u1_star,
                               sigma_sq = param_out[k,]$u2_star),
                      param = list(m_n = m_n, w_n = w_n,
                                   r_n = r_n, s_n = s_n))
    
    
    # 1st param calculation
    p1 = -1/lambda[1] * exp(-lambda[1] * (param_out[k,]$u1_ub - param_out[k,]$u1_lb))
    
    # 2nd param calculation
    p2 = -1/lambda[2] * exp(-lambda[2] * (param_out[k,]$u2_ub - param_out[k,]$u2_lb))
    
    zhat[k] = c_k[k] * p1 * p2
}

# (3) function to evaluate c_k (psi_k) -- already implemented via psi_true()


# (4) test function that extracts median from each partition




# (5) form the approximation
















