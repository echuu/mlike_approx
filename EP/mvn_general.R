

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions
## 

source("C:/Users/ericc/mlike_approx/bayes_regression/bayesLinRegHelper.R") 
# source("C:/Users/ericc/mlike_approx/paper_simulations/table2/mvn_estimators.R") 

# setwd("/mlike_approx/algo")
# source("setup.R")           # setup global environment, load in algo functions
# source("/mlike_approx/bayes_regression/bayesLinRegHelper.R") 
# source("/mlike_approx/paper_simulations/table2/mvn_estimators.R")


D = c(4) # test for smalller dimensions for now
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
J         = 2000         # number of MC samples per approximation
# ------------------------------------------------------------------------------

## sample from posterior -------------------------------------------------------

# true log marginal likelihood
# lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

# library(microbenchmark)
# microbenchmark("hml" = hml_const(1, D, u_df, J, prior))
# start_time <- Sys.time()
hml_approx = hml_const(1, D, u_df, J, prior)
# end_time <- Sys.time()

# end_time - start_time
hml_approx$const_vec       # -272.1245

(LIL = lil(prior, post))   # -272.1202





l1_norm = function(u, u_0) {
  sum(abs(u - u_0))
}


## (2) fit the regression tree via rpart()
u_rpart = rpart(psi_u ~ ., u_df)

## (3) process the fitted tree

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support) 


# param_out = u_star_cand(u_rpart, u_df, u_partition, D) # partition.R
param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R
# opt_part = param_out$optimal_part

# ----------------------------------------------------------------------
n_partitions = nrow(u_partition) # number of partitions 

# ----------------------------------------------------------------------

psi_partition = param_out %>% 
  dplyr::select(-c('leaf_id', 'psi_choice', 'logQ_cstar', 'n_obs'))

bounds = psi_partition %>% dplyr::select(-c("psi_star"))

# ------------------------------------------------------------------------------

#### extension starts here -----------------------------------------------------


### (1) find global mean
u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean


### (2) find point in each partition closest to global mean (for now)

# u_k for each partition
u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>% 
  group_by(leaf_id) %>% filter(l1_cost == min(l1_cost)) %>% 
  data.frame


K = nrow(bounds)
log_terms = numeric(K) # store terms so that we can use log-sum-exp()
G_k = numeric(K)       # store terms coming from gaussian integral
for (k in 1:K) {
  
  u_k = unname(unlist(psi_df[k,1:D]))
  # diff_k = u_k - m_k
  
  H_k = pracma::hessian(psi, u_k, prior = prior)
  H_k_inv = solve(H_k)
  lambda_k = pracma::grad(psi, u_k, prior = prior)
  b_k = H_k %*% u_k - lambda_k
  m_k = H_k_inv %*% b_k
  
  lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
  ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
  G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)

  log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
    psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
  
}

log_sum_exp(log_terms)
(LIL = lil(prior, post)) 











