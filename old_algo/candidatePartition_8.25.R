



### things to return/keep track of
### (1) partition input for the next stage of resampling
### (2) candidate_df with with leaf_id and 5 candidates
### (3) exp terms using each of the 5 candidates
### (4) log ML approximation for each of the 5 candidates
next_stage = function(curr_part, D, n_samps, params, psi_candidates) {
    
    part_id = curr_part$leaf_id        # extract partition IDs
    K = length(part_id)                # number of partitions
    
    candidate_psi = vector("list", K) # store the sub-partitions after resample
    opt_psi_part  = vector("list", K) # store the sub-partitions after resample
    u_df_k_list   = vector("list", K) # store resampled points w/ psi, psi_star
    for (k in 1:K) {

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
        # names(psi_tilde) = c("leaf_id", 
        #                      paste("psi_tilde_", 
        #                            c(1:(ncol(psi_tilde)-1)), sep = ''))
        names(psi_tilde) = c("leaf_id", 
                             paste("psi_2_", c(1:(ncol(psi_tilde)-1)), sep = ''))
        resid_partition = psi_star %>% merge(resid_partition, by = 'leaf_id')
        
        
        # ----------------------------------------------------------------------
        psi_cand_k = psi_candidates %>% dplyr::filter(leaf_id == part_id[k]) %>% 
            dplyr::select(-c('leaf_id'))
        u_df_k_cand = u_df_k %>% dplyr::mutate(leaf_id = resid_tree$where)
        # print(psi_cand_k)
        psi_tilde = data.frame(psi_tilde, psi_cand_k)
        
        u_df_k_cand = merge(u_df_k_cand, psi_tilde, by = 'leaf_id')
        # u_df_k_cand %>% head
        
        
        # sample N_k rows from the candidate data frame
        u_df_k_cand_sub = u_df_k_cand[sample(N_k_p, curr_part$n_obs[k]),]
        
        
        # ----------------------------------------------------------------------
        
        #### store (1) candidates, (2) optimal partition
        candidate_psi[[k]] = psi_tilde
        opt_psi_part[[k]]  = resid_partition
        u_df_k_list[[k]]   = u_df_k_cand_sub %>% dplyr::select(-c('leaf_id')) 
        
    } # end of for loop iterating over partitions
    
    ## all candidates
    candidate_df = do.call(rbind, candidate_psi)
    
    # print(candidate_df %>% head)
    
    ## all optimal
    optimal_df = do.call(rbind, opt_psi_part)
    u_df_resample = do.call(rbind, u_df_k_list)
    
    ## re-index the leaf IDs
    K = nrow(candidate_df)
    candidate_df = candidate_df %>% dplyr::mutate(leaf_id = 1:K)
    optimal_df = optimal_df %>% dplyr::mutate(leaf_id = 1:K)
    
    # print(optimal_df   %>%  head)
    # print(candidate_df %>%  head)
    
    return(list(candidate_psi = candidate_df,
                optimal_part  = optimal_df,
                u_df_resample = u_df_resample))
    
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
    # print(logml_approx)
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
    names(psi_tilde_df) = names(candidate_df)[-1]
    for (i in 1:n_cand) {
        # identify the psi candidate value from candidate_df
        psi_cand_i = candidate_df[,i+1]
        psi_tilde_df[,i] = psi_cand_i[part_id]
        if (stage > 1) {
            psi_tilde_df[use_prev_id,i] = psi_star[use_prev_id, stage-1]
        }
    }
    
    
    # error = apply(psi_tilde_df, 2, MAE, u_df$psi_u)
    # error = apply(psi_tilde_df, 2, logQ, u_df$psi_u)
    # print(error)
    # approx = compute_approx(part_fit)
    
    psi_star[,stage] = optimal_df$psi_star[part_id]
    if (stage > 1) {
        psi_star[,stage][use_prev_id] = psi_star[use_prev_id, stage-1]
    }
    
    # approx_error = data.frame(approx = approx, error = error)
    
    
    # return(list(approx_error = approx_error, psi_star = psi_star, 
    #             psi_tilde_df = psi_tilde_df))
    return(list(psi_star = psi_star))
    
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
logml_old = function(D, u_df, J, param) {
    
    n_stage = 2
    n_samps = 10
    params = param
    approx = vector("list", n_stage)
    psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it
    
    logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
    prev_part = logml_approx$param_out$optimal_part # first stage partition
    psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)
    
    
    
    # --------------------------------------------------------------------------
    
    n_cand = 5
    psi_candidates = logml_approx$param_out$candidate_psi
    names(psi_candidates) = c('leaf_id', paste("psi_1", 1:n_cand, sep = '_'))
    
    ## create psi_star column -- using join() preserves original ordering 
    ## of u_df, which is important when later constructing the psi_tilde 
    ## dataframe -> psi_star is pulled from this dataframe
    u_df_star = u_df %>% 
        dplyr::mutate(leaf_id = logml_approx$u_rpart$where) %>% 
        plyr::join(prev_part %>% dplyr::select(c('leaf_id', 'psi_star')),
                   by = 'leaf_id') %>% 
        dplyr::select(-c('leaf_id'))
    
    
    # --------------------------------------------------------------------------
    
    # approx[[1]] = psi_update$approx_error
    psi_star = psi_update$psi_star
    
    for (s in 2:n_stage) {
        stage_s_part = next_stage(prev_part, D, n_samps, params, psi_candidates)
        psi_update = compute_weights(u_df, stage_s_part, s, psi_star)
        
        
        approx[[s]] = psi_update$approx_error
        
        prev_part = stage_s_part$optimal_part
        psi_star = psi_update$psi_star # used in next iteration
        n_samps = n_samps / 2
    }
    
    # print(do.call(rbind, approx))
    all_approx = do.call(rbind, approx) %>% 
        data.frame(stage = sort(rep(1:n_stage, n_cand), decreasing = TRUE))
    all_approx = all_approx %>% dplyr::mutate(wt = exp(-error/(2*D))) %>% 
        dplyr::mutate(norm_wt = wt / sum(wt))
    all_approx = all_approx %>% dplyr::mutate(wt2 = exp(-error/(2*sqrt(D)))) %>% 
        dplyr::mutate(norm_wt2 = wt2 / sum(wt2))
    
    wt_approx1 = with(all_approx, sum(approx * norm_wt))
    wt_approx2 = with(all_approx, sum(approx * norm_wt2))
    avg_approx = mean(all_approx$approx)
    
    return(list(wt_approx1 = wt_approx1, wt_approx2 = wt_approx2,
                avg_approx = avg_approx,
                all_approx = all_approx,
                psi_tilde_df = psi_update$psi_tilde_df,
                u_df_ss = stage_s_part$u_df_resample))
}






logml_call = function(D, u_df, J, param, n_stage = 2, n_samps = 10) {
    
    params = param
    approx = vector("list", n_stage)
    psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it
    
    logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
    prev_part = logml_approx$param_out$optimal_part # first stage partition
    # psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)
    
    # --------------------------------------------------------------------------
    
    n_cand = 5
    psi_candidates = logml_approx$param_out$candidate_psi
    names(psi_candidates) = c('leaf_id', paste("psi_1", 1:n_cand, sep = '_'))
    
    ## create psi_star column -- using join() preserves original ordering 
    ## of u_df, which is important when later constructing the psi_tilde 
    ## dataframe -> psi_star is pulled from this dataframe
    u_df_star = u_df %>% 
        dplyr::mutate(leaf_id = logml_approx$u_rpart$where) %>% 
        plyr::join(prev_part %>% dplyr::select(c('leaf_id', 'psi_star')),
                   by = 'leaf_id') %>% 
        dplyr::select(-c('leaf_id'))
    # psi_star = psi_update$psi_star
    stage_s_part = next_stage(prev_part, D, n_samps, params, psi_candidates)
    # psi_update = compute_weights(u_df, stage_s_part, s, psi_star)
    
    optimal_df = stage_s_part$optimal_part
    candidate_df = stage_s_part$candidate_psi
    
    u_sub = u_df %>% dplyr::select(-psi_u)
    drops   = c('leaf_id', 'psi_star', 'n_obs')
    
    part = optimal_df[,!(names(optimal_df) %in% drops)]
    part_list = lapply(split(part, seq(NROW(part))), matrix_part)
    
    # for each posterior sample, determine which 2nd stage partition 
    # it belongs in (note these are indexed from 1:n_partitions)
    part_id = apply(u_sub, 1, query_partition, part_list = part_list)
    
    # for posterior samples that fall outside of the re-sampled partitions,
    # we will use the psi_star value obtained during the previous stage
    use_prev_id = which(is.na(part_id))
    
    n_cand = ncol(candidate_df) - 1 # subtract off the leaf_id column
    # each row corresponds to the same row in u_df / u_df_star / u_df_cand
    # we populate psi value for each of these rows, for each of the candidates
    # if the point u does not fall into one of the second stage partitions, i.e., 
    # it does not have a corresponding second stage psi_tilde value, we use the
    # OPTIMAL PSI VALUE chosen from the previous stage. This optimal psi value is
    # stored in u_df_star created above, where the order of u_df is preserved!
    psi_tilde_df = data.frame(matrix(0, nrow(u_df), n_cand))
    names(psi_tilde_df) = names(candidate_df)[-1]
    for (i in 1:n_cand) {
        # identify the psi candidate value from candidate_df
        # This is a vector that has length equal to the number of partitions created
        # as a result of the second stage. This vector will be repeatedly re-indexed
        # depending on which (second stage) partition u falls in
        psi_cand_i = candidate_df[,i+1]
        
        # populate the column using the corresponding candidate psi value
        psi_tilde_df[,i] = psi_cand_i[part_id]
        
        # for those points that fall outside the 2nd stage we use previous stage's
        # optimal psi value for that partition
        # psi_tilde_df[use_prev_id,i] = psi_star[use_prev_id, stage-1]
        psi_tilde_df[use_prev_id,i] = u_df_star$psi_star[use_prev_id]
    }
    
    # error = apply(psi_tilde_df, 2, logQ, u_df$psi_u) %>% 
    #     melt(value.name = 'error')
    
    log_vol = log_volume(optimal_df, drops, D)
    exp_terms_mat = -candidate_df[,-1] + log_vol
    all_approx = reshape2::melt((apply(exp_terms_mat, 2, log_sum_exp)),
                                value.name = 'approx')
    
    # all_approx
    
    # approx_fs_error = merge(all_approx, error, by = 0)
    # 
    # approx_fs_error = approx_fs_error %>% dplyr::mutate(error = -log(J) + error)
    # 
    # approx_fs_error = approx_fs_error %>% 
    #     dplyr::mutate(wt_2d = exp(-error/(2*D))) %>% 
    #     dplyr::mutate(norm_wt_2d = wt_2d / sum(wt_2d)) %>% 
    #     dplyr::select(-c('wt_2d'))
    # 
    # approx_fs_error = approx_fs_error %>% 
    #     dplyr::mutate(wt_sqrt = exp(-error/(D^2))) %>% 
    #     dplyr::mutate(norm_wt_sqrt = wt_sqrt / sum(wt_sqrt)) %>% 
    #     dplyr::select(-c('wt_sqrt'))
    
    # approx_fs_error
    # (wt_approx1 = with(approx_fs_error, sum(approx * norm_wt_2d)))
    # (wt_approx2 = with(approx_fs_error, sum(approx * norm_wt_sqrt)))
    # wt_approx1 = with(approx_fs_error, sum(approx * norm_wt_2d))
    # wt_approx2 = with(approx_fs_error, sum(approx * norm_wt_sqrt))
    
    
    # --------------------------------------------------------------------------
    
    ### combine the 1st/2nd-stage candidates together so we can compute the
    ### error at the same time
    
    num_ss = nrow(stage_s_part$u_df_resample)
    psi_ss = stage_s_part$u_df_resample[,-c(1:D)]
    
    all_df = rbind(psi_tilde_df, psi_ss %>% dplyr::select(-psi_u))
    all_psi = c(u_df$psi_u, psi_ss$psi_u)
    
    # --------------------------------------------------------------------------
    
    # print(J + num_ss)
    
    # all_error = (-log(J + num_ss) + apply(all_df, 2, logJ, all_psi)) %>% 
    #     melt(value.name = 'error')
    
    # all_error = (-log(J + num_ss) + (apply(all_df, 2, logJ, all_psi) + 
    #                                  apply(all_df, 2, logQ, all_psi))/2) %>% 
    #             melt(value.name = 'error')
    
    all_error = rbind(apply(all_df[,1:5], 2, logQ, all_psi) %>% 
                          melt(value.name = 'error'), 
                      apply(all_df[,6:10], 2, logJ, all_psi) %>%  
                          melt(value.name = 'error')) -log(J + num_ss)
    
    approx_all_error =  merge(all_approx, all_error, by = 0)
    
    approx_all_error = approx_all_error %>% 
        dplyr::mutate(wt_1 = exp(-error/(7*D))) %>% 
        dplyr::mutate(norm_wt_1 = wt_1 / sum(wt_1)) %>% 
        dplyr::select(-c('wt_1'))
    
    approx_all_error = approx_all_error %>% 
        dplyr::mutate(wt_2 = exp(-error/(8*D))) %>% 
        dplyr::mutate(norm_wt_2 = wt_2 / sum(wt_2)) %>% 
        dplyr::select(-c('wt_2'))
    
    approx_all_error = approx_all_error %>% 
        dplyr::mutate(wt_3 = exp(-error/(10*D))) %>% 
        dplyr::mutate(norm_wt_3 = wt_3 / sum(wt_3)) %>% 
        dplyr::select(-c('wt_3'))
    
    # approx_all_error
    (wt_approx1 = with(approx_all_error, sum(approx * norm_wt_1)))
    (wt_approx2 = with(approx_all_error, sum(approx * norm_wt_2)))
    (wt_approx3 = with(approx_all_error, sum(approx * norm_wt_3)))
    
    
    return(list(all_approx = approx_all_error,
                wt_approx1 = wt_approx1,
                wt_approx2 = wt_approx2,
                wt_approx3 = wt_approx3))
}


logml <- function(D, u_df, J, prior) {
    out <- tryCatch(
        {
            
            logml_call(D, u_df, J, prior)
            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            # message("Here's the original warning message:")
            # message(cond)
            return(NULL)
        },
        finally={
        }
    )    
    return(out)
}




# ### weighted average of 3 stage partitioning
# logml = function(D, u_df, J, param) {
#     
#     n_stage = 2
#     n_samps = 10
#     params = param
#     approx = vector("list", n_stage)
#     psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it
#     
#     logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
#     prev_part = logml_approx$param_out$optimal_part # first stage partition
#     
#     psi_candidates = logml_approx$param_out$candidate_psi
#     names(psi_candidates) = c('leaf_id', paste("psi_1", 1:n_cand, sep = '_'))
#     
#     
#     psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)
#     
#     approx[[1]] = psi_update$approx_error
#     psi_star = psi_update$psi_star
#     
#     for (s in 2:n_stage) {
#         stage_s_part = next_stage(prev_part, D, n_samps, params)
#         psi_update = compute_weights(u_df, stage_s_part, s, psi_star)
#         
#         approx[[s]] = psi_update$approx_error
#         
#         prev_part = stage_s_part$optimal_part
#         psi_star = psi_update$psi_star # used in next iteration
#         n_samps = n_samps / 2
#     }
#     
#     all_approx = do.call(rbind, approx) %>% 
#         data.frame(stage = sort(rep(1:n_stage, 5)))
#     all_approx = all_approx %>% dplyr::mutate(wt = exp(-error/(2*D))) %>% 
#         dplyr::mutate(norm_wt = wt / sum(wt))
#     all_approx = all_approx %>% dplyr::mutate(wt2 = exp(-error/(2))) %>% 
#         dplyr::mutate(norm_wt2 = wt2 / sum(wt2))
#     
#     wt_approx1 = with(all_approx, sum(approx * norm_wt))
#     wt_approx2 = with(all_approx, sum(approx * norm_wt2))
#     avg_approx = mean(all_approx$approx)
#     
#     return(list(wt_approx1 = wt_approx1, wt_approx2 = wt_approx2,
#                 avg_approx = avg_approx,
#                 all_approx = all_approx,
#                 psi_tilde_df = psi_update$psi_tilde_df))
# }
# 
# 
# 














