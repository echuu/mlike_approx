



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


D = c(15) # test for smalller dimensions for now
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


## compute posterior parameters ------------------------------------------------

Q_beta = 1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv = solve(Q_beta)
b = 1 / sigmasq * t(X) %*% y
mu_beta = Q_beta_inv %*% b

prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)


## algorithm settings ----------------------------------------------------------

J         = 2000          # number of MC samples per approximation
N_approx  = 1            # number of approximations to compute using algorithm
K_sims    = 1            # number of simulations to run

# ------------------------------------------------------------------------------


## sample from posterior -------------------------------------------------------

# true log marginal likelihood
lil(prior, post) 

u_samps = rmvnorm(J, mean = c(mu_beta), sigma = Q_beta_inv) %>% data.frame 
u_df = preprocess(u_samps, D, prior) # J x (D + 1) -- stored row-wise 

# hml_approx = hml(N_approx, D, u_df, J, prior) 

hml_approx = hml_const(1, D, u_df, J, prior)

hml_approx$const_vec  # -457.0106


# ------------------------------------------------------------------------------

#### model output diagnostics

hml_approx$n_partitions

part_info = hml_approx$param_out
u_df_info = hml_approx$u_df_fit

# (1.1) for each point: fitted value, leaf_id, residual

## notation: 
# psi_u     : true value, psi(u)
# psi_star  : psi(u_star), where u_star is the partition's representative point
# psi_resid : psi_u - psi_star
# psi_hat   : fitted value for psi(u) on a given partition
##

u_df_info %>% head

# (1.2) for each partition: 'median' points, fitted value, upper/lower bounds,
# number of observations in partition
part_info 


# compare the psi values from tree vs. psi values evaluated at 'representative'
# point of each partition
part_psi = part_info %>% 
    dplyr::select(leaf_id, psi_hat) %>% 
    merge(u_df_info %>% dplyr::select(leaf_id, psi_star) %>% unique, 
          by = 'leaf_id')

part_psi

# (2) fitted value for each partition (const_approx)
# (2.1) look at the MSE for each of the partitions to see if there's one
# parition where there is significantly larger discrepancy
# (2.2) consider the number of observations in each of the partitions in
# conjunction with the MSE for each partition

library(MLmetrics)
part_mse = u_df_info %>% 
    dplyr::group_by(leaf_id) %>% 
    summarise(psi_star_mse = MSE(psi_u, psi_star)) %>% 
    merge(part_info %>% dplyr::select(leaf_id, n_obs), by = 'leaf_id')


# for each partition, display the fitted value for psi (from tree), 
# psi evaluated at the representative point, the mse associated with the
# representative point, the number of observations for that partition
psi_df = merge(part_psi, part_mse, by = 'leaf_id')

psi_df %>% arrange(psi_star_mse)

u_df_info %>% filter(leaf_id == 28)

psi_df %>% arrange(n_obs)

# compute mse for the tree's fitted values
psi_rpart = hml_approx$u_rpart
tree_dev = rpart_mse(psi_rpart, u_df)
tree_mse = tree_dev %>% group_by(leaf_id) %>% 
    summarise(psi_hat_mse = mean(dev_sq))

# gather all psi approximations (tree and algo) with associated MSEs
# into one table
psi_mse = psi_df %>% merge(tree_mse, by = 'leaf_id') %>% 
    dplyr::select(leaf_id, psi_hat, psi_hat_mse, psi_star, psi_star_mse, n_obs)


psi_mse

hml_approx$const_vec 

lil(prior, post) 







































# (3) compute approximate integral over each partition
hml_approx$const_approx






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









