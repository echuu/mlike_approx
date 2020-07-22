



# # DELL_PATH = "C:/Users/chuu/mlike_approx"
# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# # path for lenovo
# setwd(LEN_PATH)
# 
# # path for dell
# # setwd(DELL_PATH)
# 
# library(mvtnorm)
# 
# source("partition/partition.R")
# source("extractPartition.R")
# source("hybrid_approx.R")
# 


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

# library(microbenchmark)
# microbenchmark("hml" = hml_const(1, D, u_df, J, prior))

hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)

hml_approx$const_vec # -272.1245

lil(prior, post)     # -272.1202





#### D = 2 (J = 500)

## sample run (1)
# always taking max   : -224.5714
# using loss function : -223.7904
# truth               : -223.7641

## sample run (2)
# always taking max   : -217.1747
# using loss function : -216.5745
# truth               : -216.599


#### D = 4 (J = 500)

## sample run (1)
# always taking max   : -240.0592
# using loss function : -239.0824
# truth               : -238.9488

## test logMSE() function


### TODO: put the stuff below in a separate function
# can be used as a diagnostics function - using the MSE doesn't seem like it 
# will yield good results empirically. for higher dimensions, it still picks
# mean/median
# for lower dimension, the two agree on a few, but still not most of them



param_support = extractSupport(u_df, D)
rpart_obj = hml_approx$u_rpart
partition = extractPartition(hml_approx$u_rpart, param_support) 

u_df = u_df %>% mutate(leaf_id = rpart_obj$where)

psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

part_obs_tbl = table(rpart_obj$where) %>% data.frame
names(part_obs_tbl) = c("leaf_id", "n_obs")

n_partitions = length(partition_id)

# compute max for each of the partitions
psi_all = u_df %>% dplyr::group_by(leaf_id) %>% 
    summarise(psi_max  = max(psi_u), 
              psi_med  = median(psi_u), 
              psi_mean = mean(psi_u),
              psi_85   = quantile(psi_u, 0.85),
              psi_90   = quantile(psi_u, 0.90),
              psi_95   = quantile(psi_u, 0.95)) %>% 
    merge(psi_hat_df, by = 'leaf_id')

psi_long_logQ = melt(psi_all, id.vars = c("leaf_id"), 
                     value.name = "psi_star_Q", 
                     variable.name = "psi_choice_Q")

psi_long_logMSE = melt(psi_all, id.vars = c("leaf_id"), 
                       value.name = "psi_star_mse", 
                       variable.name = "psi_choice_mse")

psi_logQ_df = psi_long_logQ %>% dplyr::mutate(logQ_cstar = 0)
psi_logMSE_df = psi_long_logMSE %>% dplyr::mutate(lmse_cstar = 0)

for (k in 1:n_partitions) {
    # extract psi_u for the k-th partition
    c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u
    
    psi_logQ_df = psi_logQ_df %>% 
        mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                   sapply(psi_star_Q, logQ, c_k = c_k),
                                   logQ_cstar))
    psi_logMSE_df = psi_logMSE_df %>% 
        mutate(lmse_cstar = ifelse(leaf_id == partition_id[k],
                                   sapply(psi_star_mse, logMSE, c_k = c_k),
                                   lmse_cstar))
    
    
} # end of loop extracting representative points

# for each partition (leaf_id), subset out rows for which log(Q(c)) is min
psi_min_logQ = psi_logQ_df %>% 
    group_by(leaf_id) %>% 
    slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
    data.frame()

psi_min_logQ

psi_min_logMSE = psi_logMSE_df %>% 
    group_by(leaf_id) %>% 
    slice(which.min(lmse_cstar)) %>%  # extract rows that minimize log(Q(c))
    data.frame()

psi_min_logMSE

psi_merge = merge(psi_min_logQ, psi_min_logMSE, by = "leaf_id")

psi_merge


hml_approx$const_vec # -249.8994

lil(prior, post) # -248.8395



# ------------------------------------------------------------------------------



## TODO: consider 2-d case, 

## TODO: try the case when (# of params / # samples) = 1/3
# D / N = 1 / 3
# D / N = 1 / 4







k = 1

c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u

c_star_k = psi_all_df[psi_all_df$leaf_id == partition_id[k],]$psi_star[1]

logMSE(c_star_k, c_k)











c_k = c(1:10)
c_star_k = 10
log(1/length(c_k) * sum((exp(-c_k) - exp(-c_star_k))^2))

log(mean((exp(-c_k) - exp(-c_star_k))^2))
logMSE(c_star_k, c_k)



log((exp(-c_k[1]) - exp(-c_star_k))^2)
2 * (-c_k[1] + log1mexp(c_star_k - c_k[1]))






# ------------------------------------------------------------------------------

param_support = extractSupport(u_df, D)
rpart_obj = hml_approx$u_rpart
partition = extractPartition(hml_approx$u_rpart, param_support) 


#### inside u_star() function ####
u_df = u_df %>% mutate(leaf_id = rpart_obj$where)

psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

part_obs_tbl = table(rpart_obj$where) %>% data.frame
names(part_obs_tbl) = c("leaf_id", "n_obs")

n_partitions = length(partition_id)

# current function goes to here ------------------------------------------------

## start new u_star() code

# compute max for each of the partitions
psi_center = u_df %>% dplyr::group_by(leaf_id) %>% 
    summarise(psi_med  = median(psi_u), 
              psi_mean = mean(psi_u)) %>% 
    merge(psi_hat_df, by = 'leaf_id')

psi_quant = u_df %>% dplyr::group_by(leaf_id) %>% 
    do(data.frame(t(quantile(.$psi_u, probs = seq(0.75, 1, 0.01)))))

names(psi_quant) = c("leaf_id", paste('psi_', seq(75, 100, 1), sep = ''))

psi_all = merge(psi_center, psi_quant, by = 'leaf_id') 


psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                variable.name = "psi_choice")

psi_all_df = psi_long %>% 
    dplyr::mutate(logQ_cstar = 0)


for (k in 1:n_partitions) {
    
    
    # take psi_star to be the one that gives minimum log Q
    
    # extract psi_u for the k-th partition
    c_k = u_df[u_df$leaf_id == partition_id[k],]$psi_u
    
    psi_all_df = psi_all_df %>% 
        mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                   sapply(psi_star, logQ, c_k = c_k),
                                   logQ_cstar))
    
} # end of loop extracting representative points

# for each partition (leaf_id), subset out rows for which logQ_cstar is min

psi_df = psi_all_df %>% 
    group_by(leaf_id) %>% 
    slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
    data.frame()

## end new u_star() code


psi_df



table(psi_star$variable)



















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




# library(microbenchmark)
# microbenchmark("psi" = { a = psi(unname(unlist(u_samps[1,])), prior) },
#                "psi1" = { b = psi1(unname(unlist(u_samps[1,])), prior) })









