
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
## 

source("C:/Users/ericc/mlike_approx/bayes_regression/bayesLinRegHelper.R") 
source("C:/Users/ericc/mlike_approx/paper_simulations/table2/mvn_estimators.R") 

# setwd("/mlike_approx/algo")
# source("setup.R")           # setup global environment, load in algo functions
# source("/mlike_approx/bayes_regression/bayesLinRegHelper.R") 
# source("/mlike_approx/paper_simulations/table2/mvn_estimators.R")


D = c(20) # test for smalller dimensions for now
N = c(100) # for testing -- comment this line to perform ext. analysis


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
J = 1000         # number of MC samples per approximation
# ------------------------------------------------------------------------------

## sample from posterior -------------------------------------------------------

# true log marginal likelihood
# lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 


hme_approx(u_df, prior, J, D, N)
came_approx(u_df, prior, post, J, D)
hml_simple(D, u_df, J, prior)$zhat
hybrid_ml(D, u_df, J, prior)$zhat


(LIL = lil(prior, post))   # -272.1202



#### begin simulations ---------------------------------------------------------

B = 100 # number of replications
J = 1000 # number of MCMC samples per replication

hme     = numeric(B) # store harmonic mean estiator
hyb     = numeric(B) # store hybrid estimator
bridge  = numeric(B) # store arithmetic mean estimator
came    = numeric(B) # corrected arithmetic mean estimator

set.seed(1)
for (i in 1:B) {
    
    # sample from posterior
    u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
    u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 
    
    #### (1) hybrid estimator
    hybrid_v0 = hml_simple(D, u_df, J, prior)
    # hybrid_v0$param_out
    # hybrid_v0$zhat
    hyb[i] = hybrid_v0$zhat
    
    #### (2) harmonic mean estimator
    hme[i] = hme_approx(u_df, prior, J, D, N)
    
    #### (3) corrected arithmetic mean estimator (IS)
    came[i] = came_approx(u_df, prior, post, J, D)
    
    lb <- c(rep(-Inf, D))
    ub <- c(rep(Inf, D))
    u_samp = as.matrix(u_samps)
    colnames(u_samp) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridge_sampler(samples = u_samp,
                                   log_posterior = log_density,
                                   data = prior, lb = lb, ub = ub, silent = TRUE)
    bridge[i] = bridge_result$logml
    
    if (i %% 10 == 0) { print(paste("iter:", i)) }
}

approx = data.frame(LIL, 
                    hme = hme[hme!=0],
                    hyb = hyb[hyb!=0],
                    bridge = bridge[bridge!=0],
                    came = came[came!=0])

(error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
                    mae = colMeans(abs(LIL - approx)),
                    rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3))


approx_df = data.frame(mcmc = 1:B, hme = hme, came = came, hyb = hyb)
approx_long = melt(approx_df, id.vars = 'mcmc')

ggplot(approx_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 0.9)

### mean, error
(LIL = lil(prior, post))


round(mean(hyb), 3)
round(mean(hme), 3)
# mean(ame)
round(mean(came), 3)

round(mean(LIL - hyb), 3)
round(mean(LIL - hme), 3)
# mean(LIL - ame)
round(mean(LIL - came), 3)

round(sqrt(mean((LIL - hyb)^2)), 3)
round(sqrt(mean((LIL - hme)^2)), 3)
# sqrt(mean((LIL - ame)^2))
round(sqrt(mean((LIL - came)^2)), 3)















