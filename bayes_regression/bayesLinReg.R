



# DELL_PATH = "C:/Users/chuu/mlike_approx"
LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
setwd(LEN_PATH)

# path for dell
# setwd(DELL_PATH)

library(mvtnorm)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")



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


D = c(3) # test for smalller dimensions for now
N = c(200) # for testing -- comment this line to perform ext. analysis


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

J         = 500          # number of MC samples per approximation
N_approx  = 1            # number of approximations to compute using algorithm
K_sims    = 1            # number of simulations to run

# ------------------------------------------------------------------------------


## sample from posterior -------------------------------------------------------

# true log marginal likelihood
lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

hml_approx = hml(N_approx, D, u_df, J, prior) 

hml_approx = hml_const(1, D, u_df, J, prior)

hml_approx$const_vec  # -457.1004
hml_approx$hybrid_vec




# ------------------------------------------------------------------------------

# testing the stabilizing computation for the D-dim integral

k_part = nrow(hml_approx$lambda)

ind = 2
hml_approx$ck_3[ind]


upper = hml_approx$partition$u18_ub[k_part]
lower = hml_approx$partition$u18_lb[k_part]

l_k_d = hml_approx$lambda[k_part,ind]

# log(-1 / hml_approx$lambda[6,ind] * 
#     exp(- hml_approx$lambda[6,ind] * hml_approx$partition$u18_ub[6]) - 
#     exp(- hml_approx$lambda[6,ind] * hml_approx$partition$u18_lb[6]))

- l_k_d * upper + log(- 1 / l_k_d * (1 - exp(-l_k_d * lower + l_k_d * upper)))

hml_approx$ck_3[ind]



# library(microbenchmark)
# microbenchmark("psi" = { a = psi(unname(unlist(u_samps[1,])), prior) },
#                "psi1" = { b = psi1(unname(unlist(u_samps[1,])), prior) })









