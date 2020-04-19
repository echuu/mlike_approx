



# setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

# load this LAST to overwrite def preprocess()
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
options(mc.cores = parallel::detectCores()) 


J         = 300          # number of MC samples per approximation
N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


K_sims = 1               # num of simulations to run FOR EACH N in N_vec
D = 10
N = 100 # for testing -- comment this line to perform ext. analysis


set.seed(123)


p       = D - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
a_0     = 2 / 2            # shape param for sigmasq
b_0     = 1 / 2            # scale param 
beta    = sample(-10:10, p, replace = T)
sigmasq = 4                # true variance (1 x 1) 

I_p = diag(1, p)           # (p x p) identity matrix
I_N = diag(1, N)           # (N x N) identity matrix

X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))

y = X %*% beta + eps # (N x 1) response vector
# ------------------------------------------------------------------


## compute posterior parameters ------------------------------------
V_beta_inv = solve(V_beta)
V_star_inv = t(X) %*% X + V_beta_inv

V_star  = solve(V_star_inv)                                # (p x p)
mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
a_n =  a_0 + N / 2 
b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                          t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                          t(mu_star) %*% V_star_inv %*% mu_star))

# compute MLE estimates for mean, variance of the regression model
ybar = X %*% mu_star
sigmasq_mle = 1 / N * sum((y - ybar)^2)

# create prior, posterior objects
prior = list(V_beta = V_beta, 
             mu_beta = mu_beta, 
             a_0 = a_0, 
             b_0 = b_0,
             y = y, X = X,
             V_beta_inv = V_beta_inv)

# store posterior parameters
post  = list(V_star  =  V_star,
             mu_star =  mu_star,
             a_n     =  a_n,
             b_n     =  b_n,
             V_star_inv = V_star_inv)

## form the approximation
post_dat = list(p = p,
                a_n = a_n, b_n = b_n, 
                mu_star = c(mu_star), V_star = V_star)

mvnig_fit = stan(file   = 'C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_sampler.stan', 
                 data    = post_dat,
                 iter    = J_iter,
                 warmup  = burn_in,
                 chains  = n_chains,
                 refresh = 0) # should give us J * N_approx draws

# use special preprocess b/c we call psi_true() 
u_df = preprocess(mvnig_fit, D, post, prior)

# old version of hml()
hml_approx = hml(1, D, u_df, J, prior) 
hml_approx$const_vec # -230.9012

# refactored version
hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx$const_vec # -230.9012

lil(y, X, prior, post)


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

hml_approx$const_vec # 

lil(y, X, prior, post)

psi_mse %>% summarise(mean(psi_star))








