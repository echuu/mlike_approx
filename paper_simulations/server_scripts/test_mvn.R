
setwd("/home/grad/ericchuu/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
source("/home/grad/ericchuu/mlike_approx/bayes_regression/bayesLinRegHelper.R")
source("/home/grad/ericchuu/mlike_approx/paper_sims/mvn_estimators.R")


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
J         = 5000         # number of MC samples per approximation
# ------------------------------------------------------------------------------

## sample from posterior -------------------------------------------------------

# true log marginal likelihood
# lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

# library(microbenchmark)
# microbenchmark("hml" = hml_const(1, D, u_df, J, prior))
start_time <- Sys.time()
hml_approx = hml_const(1, D, u_df, J, prior)
end_time <- Sys.time()
end_time - start_time

hml_approx$const_vec       # -272.1245
n_samps = 10

# for the partition learned from prev fitted tree, extract the partition id and
# the optimal value of psi for this partition
og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))
ss_part = fit_resid(og_part, D, n_samps, prior)
ts_part = fit_resid(ss_part, D, n_samps / 2, prior)

log_sum_exp(unlist(compute_expterms(ss_part, D)))
log_sum_exp(unlist(compute_expterms(ts_part, D)))

hme = hme_approx(u_df, prior, J, D, N)

hml_approx$const_vec       # -272.1245
hme
came_approx(u_df, hml_approx, prior, post, J, D)
(LIL = lil(prior, post))   # -272.1202

mean(c(log_sum_exp(unlist(compute_expterms(ss_part, D))),
       log_sum_exp(unlist(compute_expterms(ts_part, D))),
       hml_approx$const_vec))





