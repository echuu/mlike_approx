
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
## 

source("C:/Users/ericc/mlike_approx/bayes_regression/bayesLinRegHelper.R") 


D = c(2) # test for smalller dimensions for now
N = c(50) # for testing -- comment this line to perform ext. analysis


## priors ----------------------------------------------------------------------
set.seed(1)
mu_0 = rep(0, D)      # prior mean for beta
tau  = 1 / 4          # precision: inverse of variance
sigmasq = 4           # true variance (1 x 1) 

## true beta -------------------------------------------------------------------

beta = sample(-10:10, D, replace = T) 

## simulate regression data ----------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector

data = list(X = X, y = y)


## compute posterior parameters ------------------------------------------------

Q_beta = 1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv = solve(Q_beta)
b = 1 / sigmasq * t(X) %*% y
mu_beta = Q_beta_inv %*% b

prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)


## algorithm settings ----------------------------------------------------------
J         = 5000         # number of MC samples per approximation
# ------------------------------------------------------------------------------

## sample from posterior -------------------------------------------------------

# true log marginal likelihood
# lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

# library(microbenchmark)
# microbenchmark("hml" = hml_const(1, D, u_df, J, prior))

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


hml_approx = hml_const(1, D, u_df, J, prior)
hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)
hme = hme_approx(u_df, prior, J, D, N)


hml_approx$const_vec       # -272.1245
hme
LIL = lil(prior, post)     # -272.1202


B = 1000 # number of replications
J = 1000 # number of MCMC samples per replication

hme  = numeric(B) # store harmonic mean estiator
hyb  = numeric(B) # store hybrid estimator
ame  = numeric(B) # store arithmetic mean estimator
came = numeric(B) # corrected arithmetic mean estimator

# other estimators: chib's, bridge, more recent ones?

for (b in 1:B) {
    
    # sample from posterior
    u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
    u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 
    
    #### (1) hybrid estimator
    hml_approx = hml_const(1, D, u_df, J, prior)
    hyb[b] = hml_approx$const_vec
    
    #### (2) harmonic mean estimator
    hme[b] = hme_approx(u_df, prior, J, D, N)
    
    #### (3) corrected arithmetic mean estimator (IS)
    
    
    if (b %% 100 == 0) { print(paste("iter:", b)) }
}



approx_df = data.frame(mcmc = 1:B, hme = hme, hyb = hyb)
approx_long = melt(approx_df, id.vars = 'mcmc')

ggplot(approx_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)

### mean, error
mean(hyb)
mean(hme)
mean(ame)
mean(came)

mean(LIL - hyb)
mean(LIL - hme)
mean(LIL - ame)
mean(LIL - came)

sqrt(mean((LIL - hyb)^2))
sqrt(mean((LIL - hme)^2))
sqrt(mean((LIL - ame)^2))
sqrt(mean((LIL - came)^2))















