
## 7/15
## generate figure for partition + re-partition

## multivariate normal case


## refactored code 4/15
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
## 


source("C:/Users/ericc/mlike_approx/bayes_regression/bayesLinRegHelper.R") 


# STAN sampler settings --------------------------------------------------------
# J         = 500          # number of MC samples per approximation
# N_approx  = 1            # number of approximations
# burn_in   = 2000         # number of burn in draws
# n_chains  = 4            # number of markov chains to run
# stan_seed = 123          # seed
# 
# J_iter = 1 / n_chains * N_approx * J + burn_in 
# K_sims = 1               # num of simulations to run FOR EACH N in N_vec


D = c(2) # test for smalller dimensions for now
N = c(400) # for testing -- comment this line to perform ext. analysis


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


## compute posterior parameters ------------------------------------------------

Q_beta = 1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv = solve(Q_beta)
b = 1 / sigmasq * t(X) %*% y
mu_beta = Q_beta_inv %*% b

prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)


## algorithm settings ----------------------------------------------------------

J         = 5000         # number of MC samples per approximation
N_approx  = 1            # number of approximations to compute using algorithm
K_sims    = 1            # number of simulations to run

# ------------------------------------------------------------------------------


## sample from posterior -------------------------------------------------------

# true log marginal likelihood
# lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

# hml_approx = hml(N_approx, D, u_df, J, prior) 


hml_approx = hml_const(1, D, u_df, J, prior)
hml_approx$const_vec # -272.1245
lil(prior, post)     # -272.1202


og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

ss_part = fit_resid(og_part, D, 5, prior)
log_sum_exp(unlist(compute_expterms(ss_part, D)))



par(mfrow = c(1,2))

plot(u_df[,1], u_df[,2], pch = 20, cex = 1, 
     col = rgb(0, 0, 0, alpha = 0.15),
     xlab = '', ylab = '', main = '')

rect(hml_approx$param_out$u1_lb, hml_approx$param_out$u2_lb,
     hml_approx$param_out$u1_ub, hml_approx$param_out$u2_ub, lwd = 4,
     border = 'black')

plot(u_df[,1], u_df[,2], pch = 20, cex = 1, 
     col = rgb(0, 0, 0, alpha = 0.05),
     xlab = '', ylab = '', main = '')

rect(ss_part$u1_lb, ss_part$u2_lb,
     ss_part$u1_ub, ss_part$u2_ub, lwd = 1.8, border = 'blue')

rect(hml_approx$param_out$u1_lb, hml_approx$param_out$u2_lb,
     hml_approx$param_out$u1_ub, hml_approx$param_out$u2_ub, lwd = 2,
     border = 'black')

