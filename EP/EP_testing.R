



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

## find point in each partition closest to global mean (for now)

u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean

u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)








# compute l1 norm between each point in u_df_part and u_0

u_test = u_df_part[1:10,]
apply(u_test[,1:D], 1, l1_norm, u_0 = u_0)

l1_norm = function(u, u_0) {
    sum(abs(u - u_0))
}
u = unlist(unname(u_df_part[1,1:D]))


# use apply() to vectorize this calculation

l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>% group_by(leaf_id) %>% filter(l1_cost == min(l1_cost))

psi_df = psi_df %>% rename(psi_star = psi_u)

# 
partition_id = sort(unique(u_rpart$where)) # row id of leaf node

part_obs_tbl = table(u_rpart$where) %>% data.frame
names(part_obs_tbl) = c("leaf_id", "n_obs")


psi_df = psi_df %>% dplyr::select(psi_star, leaf_id)

partition_star = merge(psi_df, part_obs_tbl, by = 'leaf_id')
partition = u_partition %>% dplyr::select(-c('psi_hat'))
partition_star = merge(partition_star, partition)
partition_star %>% head


psi_partition = partition_star %>% 
    dplyr::select(-c('leaf_id', 'n_obs'))

bounds = psi_partition %>% dplyr::select(-c("psi_star"))

log_vol_vec = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
    log() %>% rowSums()

zhat = (-psi_partition$psi_star + log_vol_vec) %>% log_sum_exp
zhat

(LIL = lil(prior, post)) 


lb1 = bounds[1,seq(1, 2 * D, 2)] %>% unname %>% unlist
ub1 = bounds[1,seq(2, 2 * D, 2)] %>% unname %>% unlist
epmgp::pmvn(lb1, ub1, c(mu_beta), Q_beta_inv)

p_0 = TruncatedNormal::pmvnorm(c(mu_beta), Q_beta_inv, lb1, ub1)
p_0[1]


# compute G_k for each partition
K = nrow(bounds)
G_k = numeric(K)
for (k in 1:K) {
    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
    G_k[k] = epmgp::pmvn(lb, ub, c(mu_beta), Q_beta_inv, log = F)
    
    library(TruncatedNormal)
    p_0 = TruncatedNormal::pmvnorm(c(mu_beta), Q_beta_inv, lb, ub)
    p_0[1]
    
}
G_k














