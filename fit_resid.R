


# compute R(u)
# for a given u *in a specified partition*, compute psi(u) - psi_hat(u)
# where psi_hat(u) is the optimal psi value chosen in the previous layer

## start here
n_samps = 10

# for the partition learned from prev fitted tree, extract the partition id and
# the optimal value of psi for this partition
og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

set.seed(1)
ss_part = fit_resid(og_part, D, n_samps, prior)
ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
fs_part = fit_resid(ts_part, D, n_samps / 2, prior)


# truth
(true_logml = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))

# original approx
hml_approx$const_vec
log_sum_exp(unlist(compute_expterms(ss_part, D)))
log_sum_exp(unlist(compute_expterms(ts_part, D)))
log_sum_exp(unlist(compute_expterms(fs_part, D)))


# repeat the same calculations above for G-many replications -------------------

G = 50
true_logml  = numeric(G)
orig_approx = numeric(G)
ss_approx   = numeric(G)
ts_approx   = numeric(G)

for (g in 1:G) {
    
    
    X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
    eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
    y   = X %*% beta + eps                          # (N x 1) response vector
    
    # compute posterior parameters ---------------------------------------------
    Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
    Q_beta_inv =  solve(Q_beta)
    b          =  1 / sigmasq * t(X) %*% y
    mu_beta    =  Q_beta_inv %*% b
    
    # create prior, post objects to be passed into the hml algorithm 
    prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
                 Q_beta = Q_beta, b = b)
    
    lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
        0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
        1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
    
    true_logml[g] = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1])
    
    samples = data.frame(rtmvnorm(J, c(mu_beta), Q_beta_inv, 
                                  rep(0, D), rep(Inf, D)))
    
    u_df = preprocess(samples, D, prior)
    
    hml_approx = hml_const(1, D, u_df, J, prior)
    
    og_part = hml_approx$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    
    ss_part = fit_resid(og_part, D, n_samps, prior)
    ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
    
    orig_approx[g] = hml_approx$const_vec 
    
    
    ss_approx[g] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    ts_approx[g] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    # print(paste('iter ', g, '/', G, ' : ', 
    #             round(ss_approx[g], 4), ' (err: ', 
    #             round(abs(true_logml[g] - ss_approx[g]), 4), ', avg: ', 
    #             round(mean(ss_approx[1:g]), 4), '), ',
    #             round(ts_approx[g], 4), ' (err: ', 
    #             round(abs(true_logml[g] - ts_approx[g]), 4), ', avg: ', 
    #             round(mean(ts_approx[1:g]), 4), '), ', sep = ''))
    
    print(paste('iter ', g, '/', G, ' : ', 
                round(ss_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ss_approx[g]), 4), '), ',
                round(ts_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ts_approx[g]), 4), '), ', 
                sep = ''))
    
}

mean(true_logml)  # mean of the true log ML over all partitions
mean(orig_approx) # no re-partioning
mean(ss_approx)   # second stage 
mean(ts_approx)   # third stage

abs(mean(true_logml) - mean(orig_approx))
abs(mean(true_logml) - mean(ss_approx))
abs(mean(true_logml) - mean(ts_approx))







#### end of replications -------------------------------------------------------


#### fit_resid() ---------------------------------------------------------------
# curr_part : current partition: cols for leaf_id, psi_star, n_obs | ub/lbs
# D         : dimension of parameter
# n_samps   : multiplier for num (re)samples to draw from each partition
fit_resid = function(curr_part, D, n_samps, params) {
    
    
    part_id = curr_part$leaf_id        # extract partition IDs
    K = length(part_id)                # number of partitions
 
    # exp_terms = vector("list", K)    # store terms in the exp in final approx
    sub_partitions = vector("list", K) # store the sub-partitions after resample
    
    for (k in 1:K) {
        
        # print(k)
        N_k_p = curr_part$n_obs[k] * n_samps  # num of samples to draw from A_k
        
        # set of lower/upper bounds corresponding to the k-th partition
        # (1) omit the other three variables so that we only have ub/lb columns
        # (2) transform row vector into a (D x 2) matrix of lb/ub
        part_k = curr_part %>%                
            dplyr::filter(leaf_id == part_id[k]) %>%
            dplyr::select(-c(leaf_id, psi_star, n_obs)) %>% 
            unlist() %>% matrix(ncol = 2, byrow = T)
        
        # sample unif from each lower/upper bound pair to form a D-dim vector
        # lower bound in first column, upper bound in second column
        # part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        
        # sample uniformly from the D-dim partition whose lower/upper bounds
        # are defined in part_k above
        resamp_k = Matrix_runif(N_k_p, 
                                lower = part_k[,1], 
                                upper = part_k[,2]) %>% data.frame
        
        # compute psi(u) for each of the samples
        u_df_k = preprocess(resamp_k, D, params) # N_k_p x (D_u + 1)
        
        ## NEW CODE 5/29 -------------------------------------------------------
        
        # optimal value of psi_star of the original k-th partition
        c_k_star = curr_part$psi_star[k] 
        
        # compute R(u) for each of the samples; since we're in the k-th
        # partition, we can directly subtract from psi_star[k], which we have
        # computed and stored from the previous (original) call to hml_approx()
        # this will need to be modified in the (to be)-implemented func that
        # handles fitting the residuals, rather than the func values of psi(u)
        R_df_k = u_df_k %>% dplyr::mutate(R_u = psi_u - c_k_star) %>% 
            dplyr::select(-psi_u)
        
        # fit (u, R(u)) into decision tree to obtain partition
        resid_tree = rpart(R_u ~ ., R_df_k)
        
        # extract the (data-defined) support from R_df_k
        resid_support = extractSupport(R_df_k, D)
        
        # extract partition from resid_tree
        # ntoe we omit 'psi_hat' the fitted value for the residual returned from 
        # the tree, since this value results in an underestimation of the 
        # log marginal likelihood. instead, use our own objective function to 
        # choose the optimal value of R(u) for each (sub-)partition
        resid_partition = extractPartition(resid_tree, resid_support) %>% 
            dplyr::select(-psi_hat)
        
        # number of sub-partitions of partition k
        s_k = nrow(resid_partition) 
        
        # compute opt value (chosen by tree) for each sub-partition
        # e_kj_opt : leaf_id, Ru_choice, Ru_star, logJ_star, n_obs, lb/ub
        e_kj_opt = partition_opt(resid_tree, R_df_k, resid_partition, D)
        
        # extract only the optimal value (Ru_star) and partition id (leaf_id)
        # and merge this back with the partition obtained from the fitted tree
        resid_partition = e_kj_opt %>% dplyr::select(leaf_id, Ru_star) %>% 
            merge(resid_partition, by = 'leaf_id')
        
        # compute updated value for c_k_star to be used in final approx
        # update psi_star = c_k_star (prev value) + Ru_star (approx residual)
        # resid_partition : leaf_id, psi_star, lb/ub
        # (only piece of missing information now is the # of observations in 
        # each partition)
        resid_partition = resid_partition %>% 
            dplyr::mutate(psi_star = Ru_star + c_k_star) %>% 
            dplyr::select(-Ru_star)
        
        # compute areas of A_k[j], j = 1, ... , s_k
        # exp_terms[[k]] = compute_expterms(resid_partition, D)
        # prepare for next recursive call by formatting the df to match input: 
        # leaf_id, lb/ub partitions, psi_star, n_obs
        
        # obtain number of observations in each partition
        part_tbl = table(resid_tree$where) %>% data.frame
        names(part_tbl) = c("leaf_id", "n_obs")
        
        # append column of number of observations to resid_partitions
        sub_partitions[[k]] = merge(resid_partition, part_tbl, by = 'leaf_id')
        
    } # end of for loop iterating over partitions
    
    
    # combine all sub-partitions into a dataframe containing all sub-partitions
    sub_part_df = do.call(rbind, sub_partitions)
    # redefine the partition IDs
    sub_part_df = sub_part_df %>% dplyr::mutate(leaf_id = 1:nrow(sub_part_df))
    
    # note: sub_part matches the format of the input df, curr_part, and
    # can be used to recursively/iteratively call this function
    
    return(sub_part_df)
    
} # ----------------------------------------------------------------------------





## compute_expterms() ----------------------------------------------------------
## part_info: df containing leaf_id, psi_star, lb/ub of D-intervals that make 
## up each of the K partitions
compute_expterms = function(part_info, D) {
    
    # exponential term for the k-th partition, A_k: 
    # - psi_star_k + log(volume(A_k))
    
    # remove columns that are not lower/upper bounds
    # bounds = part_info %>% dplyr::select(-c(leaf_id, psi_star))
    
    # uncomment below if we want to run this function on the final subpartition
    # after exiting the loop rather than inside the loop
    bounds = part_info %>% dplyr::select(-c(leaf_id, psi_star, n_obs))
    
    # (1) compute length of each of the D intervals for each partition
    # (2) take log of each interval difference
    # (3) add up each row for log volume of each partition
    log_vol = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
        log() %>% rowSums()
    
    return(-part_info$psi_star + log_vol) # K-dim vector
} # end of compute_expterms() function -----------------------------------------




# rpart_obj    : fitted rpart object
# df           : D-dim posterior samps stored row-wise, func eval in last col
# partition_df : partition id, lb/ub of each dimension
# D            : dimension of parameter (posterior samples)
# df = R_df_k
# rpart_obj = resid_tree
# partition_df = resid_partition

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
        summarise(Ru_min  = min(R_u))
    
    
    ## include more candidates Ru* to feed into objective function
    Ru_quant = df %>% dplyr::group_by(leaf_id) %>% 
        do(data.frame(t(quantile(.$R_u, probs = seq(0.01, 0.10, 0.01)))))
    
    names(Ru_quant) = c("leaf_id", paste('Ru_', seq(1, 10, 1), sep = ''))
    
    Ru_df = merge(Ru_df, Ru_quant, by = 'leaf_id') 
    
    # --------------------------------------------------------------------------
    
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
    
    # select Ru* that minimizes objective function
    Ru_all_df = Ru_all_df %>% 
        group_by(leaf_id) %>% 
        slice(which.min(logJ_star)) %>% 
        data.frame()
    
    ## merge with the boundary of each of the partitions
    # partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id') %>%
    #     dplyr::select(-Ru_hat)
    
    partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id')
    
    # partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id')
    
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
    
    N_k_p = part_0$n_obs[k] * n_samps  # num of (re)samples to draw from part k
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





