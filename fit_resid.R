


# compute R(u)
# for a given u *in a specified partition*, compute psi(u) - psi_hat(u)
# where psi_hat(u) is the optimal psi value chosen in the previous layer

# params that will be passed into the function eventually
n_samps = 10

# for the partition learned from prev fitted tree, extract the partition id and
# the optimal value of psi for this partition
prev_info = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

# prev_info = hml_approx$param_out %>% dplyr::select(leaf_id, psi_choice, n_obs)
part_id = prev_info$leaf_id

K = length(part_id)
k = 1

exp_terms = vector("list", K) 
update_info = vector("list", K) 

for (k in 1:K) {
    print(k)
    N_k_p = prev_info$n_obs[k] * n_samps  # num of samples to draw from part k
    part_k = prev_info %>%                # set of lower/upper bounds
        dplyr::filter(leaf_id == part_id[k]) %>%
        dplyr::select(-c(leaf_id, psi_star, n_obs))
    
    # sample uniformly from each lower/upper bound pair to form a D-dim vector
    # lower bound in first column, upper bound in second column
    part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
    
    # sample uniformly from the D-dim partition whose lower/upper bounds
    # are defined in part_k_long above
    resamp_k = Matrix_runif(N_k_p, 
                            lower = part_k_long[,1],
                            upper = part_k_long[,2]) %>% data.frame
    
    # compute psi(u) for each of the samples
    u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
    
    ## NEW CODE 5/29 -----------------------------------------------------------
    
    # optimal value of psi_star of the original k-th partition
    c_k_star = prev_info$psi_star[k] 
    
    # TODO: compute R(u) for each of the samples; since we're in the k-th
    # partition, we can directly subtract from psi_star[k], which we have
    # computed and stored from the previous (original) call to hml_approx()
    # this will need to be modified in the (to be)-implemented function that
    # handles fitting the residuals, rather than the function values of psi(u)
    R_df_k = u_df_k %>% dplyr::mutate(R_u = psi_u - c_k_star) %>% 
        dplyr::select(-psi_u)
    
    # TODO: fit (u, R(u)) into decision tree to obtain partition
    resid_tree = rpart(R_u ~ ., R_df_k)
    
    # TODO: extract the (data-defined) support from R_df_k
    resid_support = extractSupport(R_df_k, D)
    
    # TODO: extract partition from resid_tree
    resid_partition = extractPartition(resid_tree, resid_support) 
    
    # TODO: for each of the partitions, obtain the optimal value as returned
    # from the tree (can use previously defined objective function)
    # this is stored in the psi_hat column of 'resid_partition'
    
    # number of sub-partitions of partition k
    s_k = nrow(resid_partition) 
    
    # opt value (chosen by tree) for each sub-partition
    # e_kj_star = resid_partition$psi_hat # (s_k x 1)-dim vector 
    
    e_kj_star = resid_partition %>% dplyr::select(leaf_id, psi_hat)
    
    # ss_partition[[k]] = c_k_approx
    e_kj_opt = partition_opt(resid_tree, R_df_k, resid_partition, D)
    
    resid_partition = e_kj_opt %>% dplyr::select(leaf_id, Ru_star) %>% 
        merge(resid_partition, by = 'leaf_id')
    
    # compute updated value for c_k_star to be used in final appro
    # c_k_star_new = c_k_star + e_kj_star$psi_hat # TODO: replace w/ overwrite later
    
    resid_partition = resid_partition %>% dplyr::rename(psi_star = psi_hat) %>% 
        dplyr::mutate(psi_star = Ru_star + c_k_star)
    
    # TODO: compute areas of A_k[j], j = 1, ... , s_k
    # using only columns that are lb/ub, subtract lb from ub
    resid_bounds = resid_partition %>% 
        dplyr::select(-c(leaf_id, Ru_star, psi_star))
    
    part_len = resid_bounds[seq(2,ncol(resid_bounds),2)] - 
        resid_bounds[seq(1,ncol(resid_bounds),2)]
    log_len = log(part_len)
    log_vol = rowSums(log_len) # log volume of each of the partitions
    
    # TODO: store residual partitions (might need later?)
    
    exp_terms[[k]] = -resid_partition$psi_star + log_vol
    
    # in order for the next recursive routine to be called, need:
    # leaf_id, psi_star, lb/ub partitions, n_obs
    part_tbl = table(resid_tree$where) %>% data.frame
    names(part_tbl) = c("leaf_id", "n_obs")
    
    update_info[[k]] = merge(resid_partition, part_tbl, by = 'leaf_id')
}


log_sum_exp(unlist(exp_terms))

hml_approx$const_vec

(true_logml = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))


sub_part = do.call(rbind, update_info)

update_info[[1]]


fit_resid = function(u_df_R, D) {
    
}


# rpart_obj    : fitted rpart object
# df           : D-dim posterior samps stored row-wise, function eval in last col
# partition_df : partition id, lb/ub of each dimension
# D            : dimension of parameter (posterior samples)


df = R_df_k
rpart_obj = resid_tree
partition_df = resid_partition


partition_opt = function(rpart_obj, df, partition_df, D) {
    
    df = df %>% dplyr::mutate(leaf_id = rpart_obj$where)
    
    # extract fitted values from tree, merged with other psis later
    # R_hat_df = partition %>% dplyr::select(leaf_id, psi_hat)
    
    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node
    
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    n_partitions = length(partition_id)
    
    ### start updated code 4/27
    Ru_df = df %>% dplyr::group_by(leaf_id) %>% 
        summarise(psi_min  = min(R_u))
    
    Ru_long = melt(Ru_df, id.vars = c("leaf_id"), value.name = "Ru_star",
                    variable.name = "Ru_choice")
    
    Ru_all_df = Ru_long %>% 
        dplyr::mutate(logJ_star = 0)
    
    
    for (k in 1:n_partitions) {
        # extract psi_u for the k-th partition
        c_k = df[df$leaf_id == partition_id[k],]$R_u
        
        # compute log(Q(c_star)) for each candidate psi_star
        Ru_all_df = Ru_all_df %>% 
            dplyr::mutate(logJ_star = ifelse(leaf_id == partition_id[k],
                                             sapply(Ru_star, logJ, c_k = c_k),
                                             logJ_star))
        
    } # end of loop extracting representative points
    
    ## merge with the boundary of each of the partitions
    partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id') %>% 
        dplyr::select(-psi_hat)
    
    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')
    
    return(partition_star)
}









part_0 = hml_obj$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

orig_partition = hml_obj$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
    arrange(desc(perc))

K = length(part_set)

# initialize a list to store the vector containing the terms in exponential
# for each of the sub-partitions
# kth elmt is an s_k dim vector of terms that are to be exponentiated
# at the very end, all entries are unlisted and evaluated via log-sum-exp

exp_terms = vector("list", K) 
ck_star_list = vector("list", K)
ss_partition = vector("list", K) # store each of the hml_obj objects

perc_thresh = sort(orig_partition$perc, decreasing = T)

for (k in 1:K) {
    
    N_k_p = part_0$n_obs[k] * n_samps  # number of (re)samples to draw from part k
    part_k = part_0 %>%           # set of lower/upper bounds
        dplyr::filter(leaf_id == part_set[k]) %>%
        dplyr::select(-c(leaf_id, n_obs))
    
    # sample uniformly from each lower/upper bound pair to form a D-dim vector
    part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
    
    resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1],
                            upper = part_k_long[,2]) %>% data.frame
    
    u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
    
    c_k_approx = hml_const_mod(1, D, u_df_k, N_k_p, prior)
    
    ck_star_list[[k]] = c_k_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)
    
    ss_partition[[k]] = c_k_approx
    
    exp_terms[[k]] = c_k_approx$const_approx
    
}





