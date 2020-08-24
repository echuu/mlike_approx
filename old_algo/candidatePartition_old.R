



### things to return/keep track of
### (1) partition input for the next stage of resampling
### (2) candidate_df with with leaf_id and 5 candidates
### (3) exp terms using each of the 5 candidates
### (4) log ML approximation for each of the 5 candidates

next_stage = function(curr_part, D, n_samps, params) {
    
    part_id = curr_part$leaf_id        # extract partition IDs
    K = length(part_id)                # number of partitions
    
    candidate_psi = vector("list", K) # store the sub-partitions after resample
    opt_psi_part  = vector("list", K) # store the sub-partitions after resample
    
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
        
        ####  8/9 from here, we can just compute the log volume of each of the 
        ####  hypercubes and store them separately (everything will be done 
        ####  within this function)
        # part_volume = data.frame(leaf_id = resid_partition$leaf_id,
        #                          log_vol = log_volume(resid_partition, D))
        
        # number of sub-partitions of partition k
        s_k = nrow(resid_partition) 
        
        # compute opt value (chosen by tree) for each sub-partition
        # e_kj_opt : leaf_id, Ru_choice, Ru_star, logJ_star, n_obs, lb/ub
        
        # 8/9 udpate: need two data structures
        # (1) optimal psi_star that minimizes loss function
        # (2) candidate psi_values to be used in exponential weighting
        
        e_kj = partition_opt_update(resid_tree, R_df_k, resid_partition, D)
        # (1) obtain optimal: [leaf_id, psi_star]
        opt = e_kj$Ru_df_opt
        # (2) obtain candidates: [leaf_id, psi_1, ..., psi_5]
        cand = e_kj$Ru_df_cand
        
        # compute psi_star = c_k_star + e_kj_star, e_kj_star = Ru_star
        psi_star = cbind(leaf_id = opt$leaf_id, 
                         psi_star = opt$Ru_star + c_k_star,
                         n_obs = opt$n_obs)
        
        # compute psi_tilde using the candidate R(u) values
        psi_tilde = cbind(leaf_id = cand$leaf_id, cand[,-1] + c_k_star)
        names(psi_tilde) = c("leaf_id", 
                             paste("psi_tilde_", 
                                   c(1:(ncol(psi_tilde)-1)), sep = ''))
        resid_partition = psi_star %>% merge(resid_partition, by = 'leaf_id')
        
        #### store (1) candidates, (2) optimal partition
        candidate_psi[[k]] = psi_tilde
        opt_psi_part[[k]]  = resid_partition
        
    } # end of for loop iterating over partitions
    
    ## all candidates
    candidate_df = do.call(rbind, candidate_psi)
    
    ## all optimal
    optimal_df = do.call(rbind, opt_psi_part)
    
    ## re-index the leaf IDs
    K = nrow(candidate_df)
    candidate_df = candidate_df %>% dplyr::mutate(leaf_id = 1:K)
    optimal_df = optimal_df %>% dplyr::mutate(leaf_id = 1:K)
    
    return(list(candidate_psi = candidate_df,
                optimal_part  = optimal_df))
    
}


### check containment for each point
matrix_part = function(row_part) {
    row_part %>% matrix(ncol = 2, byrow = TRUE)
}

check_member = function(u, A_k) {
    all(A_k[,1] <= u) && all(A_k[,2] >= u)
}

## given a point, determine which of the k partitions it belongs in
query_partition = function(u, part_list) {
    ind = (1:length(part_list))[(sapply(part_list, check_member, u = u))]
    if (length(ind) == 0) {
        ind = NA
    }
    return(ind)
}


compute_approx = function(part_fit) {
    
    optimal_df = part_fit$optimal_part
    candidate_df = part_fit$candidate_psi
    
    # compute volume of each partition
    drops   = c('leaf_id', 'psi_star', 'n_obs')
    bounds  = optimal_df[,!(names(optimal_df) %in% drops)]
    log_vol = log_volume(optimal_df, drops, D)
    
    # form the exponential term in the approximation
    exp_terms_mat = -candidate_df[,-1] + log_vol
    
    # compute the final approximation via log-sum-exp trick
    logml_approx = reshape2::melt((apply(exp_terms_mat, 2, log_sum_exp)),
                                  value.name = 'approx')
    return(logml_approx)
} # end compute_approx() function


# psi_tilde_df is (J x 3) -- (stage - 1)-th column MUST BE POPULATED for 
# stage >= 2
compute_weights = function(u_df, part_fit, stage, psi_star) {
    
    u_sub = u_df %>% dplyr::select(-psi_u)
    
    drops   = c('leaf_id', 'psi_star', 'n_obs')
    
    optimal_df = part_fit$optimal_part
    candidate_df = part_fit$candidate_psi
    
    part = optimal_df[,!(names(optimal_df) %in% drops)]
    part_list = lapply(split(part, seq(NROW(part))), matrix_part)
    
    # for each posterior sample, determine which partition it belongs in
    part_id = apply(u_sub, 1, query_partition, part_list = part_list)
    
    # for posterior samples that fall outside of the re-sampled partitions,
    # we will use the psi_star value obtained during the previous stage
    use_prev_id = which(is.na(part_id))
    
    n_cand = ncol(candidate_df) - 1 # subtract off the leaf_id column
    psi_tilde_df = data.frame(matrix(0, nrow(u_df), n_cand))
    for (i in 1:n_cand) {
        # identify the psi candidate value from candidate_df
        psi_cand_i = candidate_df[,i+1]
        psi_tilde_df[,i] = psi_cand_i[part_id]
        if (stage > 1) {
            psi_tilde_df[use_prev_id,i] = psi_star[use_prev_id, stage-1]
        }
    }
    
    
    # error = apply(psi_tilde_df, 2, MAE, u_df$psi_u)
    error = apply(psi_tilde_df, 2, logQ, u_df$psi_u)
    approx = compute_approx(part_fit)
    
    psi_star[,stage] = optimal_df$psi_star[part_id]
    if (stage > 1) {
        psi_star[,stage][use_prev_id] = psi_star[use_prev_id, stage-1]
    }
    
    approx_error = data.frame(approx = approx, error = error)
    
    
    return(list(approx_error = approx_error, psi_star = psi_star, psi_tilde_df = psi_tilde_df))
    
}




partition_opt_update = function(rpart_obj, df, partition_df, D) {
    
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
    Ru_quant = df %>% dplyr::group_by(leaf_id) %>% 
        do(data.frame(t(quantile(.$R_u, probs = seq(0.02, 0.20, 0.06)))))
    
    names(Ru_quant) = c("leaf_id", paste('Ru_', seq(2, 20, 6), sep = ''))
    
    ## include more candidates Ru* to feed into objective function
    # Ru_quant = df %>% dplyr::group_by(leaf_id) %>% 
    #     do(data.frame(t(quantile(.$R_u, probs = seq(0.10, 0.20, 0.01)))))
    
    # names(Ru_quant) = c("leaf_id", paste('Ru_', seq(10, 20, 1), sep = ''))
    
    # candidates to be returned 
    Ru_df_cand  = merge(Ru_df, Ru_quant, by = 'leaf_id') 
    
    # find optimal R(u) for the next stage
    Ru_df = Ru_df_cand
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
    partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id')
    
    # partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id')
    # partition_star = merge(Ru_df, partition_df, by = 'leaf_id') 
    # 
    # partition_star = merge(Ru_all_df, partition_df, by = 'leaf_id')
    # 
    # # append the number of observations for each leaf node to the right
    # # this is later used to determine the type of approximation to use
    partition_star = merge(partition_star, part_obs_tbl, by = 'leaf_id')
    
    return(list(Ru_df_cand = Ru_df_cand,
                Ru_df_opt  = partition_star))
    
} # end of partition_opt() function --------------------------------------------





## compute_expterms() ----------------------------------------------------------
# part_info: df containing leaf_id, psi_star, lb/ub of D-intervals that make 
#            up each of the K partitions
# D        : dimension of parameter
log_volume = function(part_info, drops, D) {
    
    
    bounds = part_info[,!(names(part_info) %in% drops)]
    
    log_vol = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
        log() %>% rowSums()
    
    # (1) compute length of each of the D intervals for each partition
    # (2) take log of each interval difference
    # (3) add up each row for log volume of each partition
    log_vol = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
        log() %>% rowSums()
    
    return(log_vol) # K-dim vector
} # end of compute_expterms() function -----------------------------------------





### weighted average of 3 stage partitioning
logml = function(D, u_df, J, param) {
    
    n_stage = 2
    n_samps = 10
    params = param
    approx = vector("list", n_stage)
    psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it
    
    logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
    prev_part = logml_approx$param_out$optimal_part # first stage partition
    psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)
    
    approx[[1]] = psi_update$approx_error
    psi_star = psi_update$psi_star
    
    for (s in 2:n_stage) {
        stage_s_part = next_stage(prev_part, D, n_samps, params)
        psi_update = compute_weights(u_df, stage_s_part, s, psi_star)
        
        approx[[s]] = psi_update$approx_error
        
        prev_part = stage_s_part$optimal_part
        psi_star = psi_update$psi_star # used in next iteration
        n_samps = n_samps / 2
    }
    
    all_approx = do.call(rbind, approx) %>% 
        data.frame(stage = sort(rep(1:n_stage, 5)))
    all_approx = all_approx %>% dplyr::mutate(wt = exp(-error/(2*D))) %>% 
        dplyr::mutate(norm_wt = wt / sum(wt))
    all_approx = all_approx %>% dplyr::mutate(wt2 = exp(-error/(2))) %>% 
        dplyr::mutate(norm_wt2 = wt2 / sum(wt2))
    
    wt_approx1 = with(all_approx, sum(approx * norm_wt))
    wt_approx2 = with(all_approx, sum(approx * norm_wt2))
    avg_approx = mean(all_approx$approx)
    
    return(list(wt_approx1 = wt_approx1, wt_approx2 = wt_approx2,
                avg_approx = avg_approx,
                all_approx = all_approx,
                psi_tilde_df = psi_update$psi_tilde_df))
}






