

param = prior
logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
prev_part = logml_approx$param_out$optimal_part # first stage partition

u_rpart = rpart(psi_u ~ ., u_df)

## (3) process the fitted tree

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support) 

# ------------------------------------------------------------------------------
param_out = u_star_cand(u_rpart, u_df, u_partition, D) # partition.R
u_df_part = u_df %>% mutate(leaf_id = u_rpart$where)

rpart_obj = u_rpart
partition = u_partition

curr_part = param_out$optimal_part
part_id = curr_part$leaf_id        # extract partition IDs
K = length(part_id)                # number of partitions

candidate_psi = vector("list", K) # store the sub-partitions after resample
opt_psi_part  = vector("list", K) # store the sub-partitions after resample
u_df_k_list   = vector("list", K) # store resampled points w/ psi, psi_star

for (k in 1:K) {
    N_k_p = curr_part$n_obs[k]  # num of samples to draw from A_k
    
    # set of lower/upper bounds corresponding to the k-th partition
    # (1) omit the other three variables so that we only have ub/lb columns
    # (2) transform row vector into a (D x 2) matrix of lb/ub
    part_k = curr_part %>%                
        dplyr::filter(leaf_id == part_id[k]) %>%
        dplyr::select(-c(leaf_id, psi_star, n_obs)) %>% 
        unlist() %>% matrix(ncol = 2, byrow = T)
    
    # sample uniformly from the D-dim partition whose lower/upper bounds
    # are defined in part_k above
    resamp_k = Matrix_runif(N_k_p, 
                            lower = part_k[,1], 
                            upper = part_k[,2]) %>% data.frame
    
    # compute psi(u) for each of the samples
    u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
    u_df_k = u_df_k %>% dplyr::mutate(leaf_id = part_id[k])
    
    u_df_k_list[[k]] = u_df_k
}

u_df_unif = do.call(rbind, u_df_k_list)

psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)

partition_id = sort(unique(rpart_obj$where)) # row id of leaf node

part_obs_tbl = table(rpart_obj$where) %>% data.frame
names(part_obs_tbl) = c("leaf_id", "n_obs")

n_partitions = length(partition_id)

psi_center = u_df_unif %>% dplyr::group_by(leaf_id) %>% 
    summarise(psi_max  = max(psi_u))

# psi_quant = u_df %>% dplyr::group_by(leaf_id) %>% 
#     do(data.frame(t(quantile(.$psi_u, probs = seq(0.85, 1, 0.05)))))
psi_quant = u_df_unif %>% dplyr::group_by(leaf_id) %>%
    do(data.frame(t(quantile(.$psi_u, probs = seq(0.9996, 0.9999, 0.00005)))))

# psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
#     do(data.frame(t(quantile(.$psi_u, 
#                              probs = seq(0.5, 0.99, length.out = 7)))))

names(psi_quant) = c("leaf_id", paste('psi_', seq(1, 7, 1), sep = ''))

psi_all = merge(psi_center, psi_quant, by = 'leaf_id') 

psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                variable.name = "psi_choice")

psi_all_df = psi_long %>% 
    dplyr::mutate(logQ_cstar = 0)


for (k in 1:n_partitions) {
    # extract psi_u for the k-th partition
    c_k = u_df_unif[u_df_unif$leaf_id == partition_id[k],]$psi_u
    
    # compute log(Q(c_star)) for each candidate psi_star
    psi_all_df = psi_all_df %>% 
        mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                   sapply(psi_star, logJ, c_k = c_k),
                                   logQ_cstar))
    
} # end of loop extracting representative points

# for each partition (leaf_id), subset out rows for which log(Q(c)) is min
psi_df = psi_all_df %>% 
    group_by(leaf_id) %>% 
    slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
    data.frame()

## merge with the boundary of each of the partitions
partition_star = merge(psi_df, partition, by = 'leaf_id') %>% 
    dplyr::select(-psi_hat)

# append the number of observations for each leaf node to the right
# this is later used to determine the type of approximation to use
partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')

# this partition will be fed into any subsequent resampling functions;
# omits extra information pertaining to loss functions, optimal psi choice 
optimal_part = partition_star %>% 
    dplyr::select(-c('psi_choice', 'logQ_cstar'))
# ------------------------------------------------------------------------------


param_out = u_star_unif(u_rpart, u_df, u_partition, D)

rpart_obj = u_rpart
partition = u_partition
n_params = D



u_star_unif = function(rpart_obj, u_df, partition, n_params) {

    u_df_part = u_df %>% mutate(leaf_id = rpart_obj$where)
    
    # curr_part = param_out$optimal_part
    
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    curr_part = merge(partition, part_obs_tbl, by = 'leaf_id')
    
    part_id = partition$leaf_id        # extract partition IDs
    K = length(part_id)                # number of partitions
    
    candidate_psi = vector("list", K) # store the sub-partitions after resample
    opt_psi_part  = vector("list", K) # store the sub-partitions after resample
    u_df_k_list   = vector("list", K) # store resampled points w/ psi, psi_star
    
    for (k in 1:K) {
        N_k_p = curr_part$n_obs[k] * 10 # num of samples to draw from A_k
        
        # set of lower/upper bounds corresponding to the k-th partition
        # (1) omit the other three variables so that we only have ub/lb columns
        # (2) transform row vector into a (D x 2) matrix of lb/ub
        part_k = curr_part %>%                
            dplyr::filter(leaf_id == part_id[k]) %>%
            dplyr::select(-c(leaf_id, psi_hat, n_obs)) %>% 
            unlist() %>% matrix(ncol = 2, byrow = T)
        
        # sample uniformly from the D-dim partition whose lower/upper bounds
        # are defined in part_k above
        resamp_k = Matrix_runif(N_k_p, 
                                lower = part_k[,1], 
                                upper = part_k[,2]) %>% data.frame
        
        # compute psi(u) for each of the samples
        u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
        u_df_k = u_df_k %>% dplyr::mutate(leaf_id = part_id[k])
        
        u_df_k_list[[k]] = u_df_k
    }
    
    u_df_unif = do.call(rbind, u_df_k_list)
    
    psi_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)
    
    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node
    
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    n_partitions = length(partition_id)
    
    psi_center = u_df_unif %>% dplyr::group_by(leaf_id) %>% 
        summarise(psi_med  = median(psi_u))
    
    # psi_quant = u_df %>% dplyr::group_by(leaf_id) %>% 
    #     do(data.frame(t(quantile(.$psi_u, probs = seq(0.85, 1, 0.05)))))
    # psi_quant = u_df_unif %>% dplyr::group_by(leaf_id) %>%
    #     do(data.frame(t(quantile(.$psi_u, 
    #                              probs = seq(0.9996, 0.9999, 0.00005)))))
    
    psi_quant = u_df_unif %>% dplyr::group_by(leaf_id) %>%
        do(data.frame(t(quantile(.$psi_u, 
                                 probs = seq(0.01, 0.9, length.out = 7)))))
    
    # psi_quant = u_df %>% dplyr::group_by(leaf_id) %>%
    #     do(data.frame(t(quantile(.$psi_u, 
    #                              probs = seq(0.5, 0.99, length.out = 7)))))
    
    names(psi_quant) = c("leaf_id", paste('psi_', seq(1, 7, 1), sep = ''))
    
    psi_all = merge(psi_center, psi_quant, by = 'leaf_id') 
    
    psi_long = melt(psi_all, id.vars = c("leaf_id"), value.name = "psi_star",
                    variable.name = "psi_choice")
    
    psi_all_df = psi_long %>% 
        dplyr::mutate(logQ_cstar = 0)
    
    
    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = u_df_unif[u_df_unif$leaf_id == partition_id[k],]$psi_u
        
        # compute log(Q(c_star)) for each candidate psi_star
        psi_all_df = psi_all_df %>% 
            mutate(logQ_cstar = ifelse(leaf_id == partition_id[k],
                                       sapply(psi_star, logJ, c_k = c_k),
                                       logQ_cstar))
        
    } # end of loop extracting representative points
    
    # for each partition (leaf_id), subset out rows for which log(Q(c)) is min
    psi_df = psi_all_df %>% 
        group_by(leaf_id) %>% 
        slice(which.min(logQ_cstar)) %>%  # extract rows that minimize log(Q(c))
        data.frame()
    
    ## merge with the boundary of each of the partitions
    partition_star = merge(psi_df, partition, by = 'leaf_id') %>% 
        dplyr::select(-psi_hat)
    
    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')
    
    # this partition will be fed into any subsequent resampling functions;
    # omits extra information pertaining to loss functions, optimal psi choice 
    optimal_part = partition_star %>% 
        dplyr::select(-c('psi_choice', 'logQ_cstar'))
    
    return(list(partition_star = partition_star,
                candidate_psi  = psi_all,
                optimal_part   = optimal_part))
    
    
}

# ------------------------------------------------------------------------------

## run the code in trunc_final.R


library(TruncatedNormal)
library(tmg)
library(mvtnorm)
library(Rcpp)

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")
sourceCpp("C:/Users/ericc/mlike_approx/speedup/trunc_psi.cpp")


param = prior
# sample from posterior
samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
u_df = preprocess(data.frame(samples), D, prior)

hml_approx = hml_const(1, D, u_df, J, prior)
hml_approx$const_vec
# 
# logml_approx = hml_approx(D, u_df, J, param)     # fit first stage
# prev_part = logml_approx$param_out$optimal_part # first stage partition

u_rpart = rpart(psi_u ~ ., u_df)

## (3) process the fitted tree

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support) 


plot(u_df[,1], u_df[,2], pch = 20, cex = 1, 
     col = rgb(0, 0, 0, alpha = 0.08),
     xlab = '', ylab = '', main = '')

rect(hml_approx$param_out$u1_lb, hml_approx$param_out$u2_lb,
     hml_approx$param_out$u1_ub, hml_approx$param_out$u2_ub, lwd = 2,
     border = 'black')

part = hml_approx$param_out %>% 
    dplyr::select('u1_lb':paste("u", D, "_ub", sep = ''))


# for each partition set, construct the convex hull 

hml_approx$u_df_fit %>% head

# leaf id, u_1, ... , u_D, psi(u)
u_df_part = hml_approx$u_df_fit %>% dplyr::select(c('leaf_id':'psi_u'))

table(u_df_part$leaf_id)

# for each leaf node, construct a convex hull
u_df_pts_k = u_df_part %>% dplyr::filter(leaf_id == 21) %>% 
    dplyr::select(-c('psi_u', 'leaf_id')) %>% 
    as.matrix

library(geometry)
library(uniformly)
u_df_k_tri = delaunayn(u_df_pts_k)
u_df_k_resample = runif_in_simplex(nrow(u_df_pts_k), u_df_k_tri)


d <- c(-1,1)
pc <- as.matrix(rbind(expand.grid(d,d,d),0))
tc <- delaunayn(pc)

install.packages("alphahull")
library("alphahull")

alpha = 0.2
ahull_obj = ahull(u_df_pts_k, alpha = alpha)







