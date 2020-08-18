


# og_part = hml_approx$param_out %>%
#     dplyr::select(-c(psi_choice, logQ_cstar))
# ss_part = fit_resid(og_part, D, n_samps, prior)
# 
# ss_part %>% head
# 
# log_sum_exp(unlist(compute_expterms(ss_part, D)))
# 
# log_sum_exp(unlist(compute_expterms(sub_part_df, D)))


# ------------------------------------------------------------------------------

### updated version of partition_opt() function that gives additional candidates
### for optimal value over each partition

# rpart_obj = resid_tree
# df = R_df_k
# partition_df = resid_partition





# ------------------------------------------------------------------------------

# 
# params = prior
# curr_part = og_part
# part_id = curr_part$leaf_id        # extract partition IDs
# K = length(part_id)                # number of partitions
# 
# 
# ### third stage stuff fuck this is a mess --------------------------------------
# 
# 
# curr_part = part2_opt
# part_id = curr_part$leaf_id        # extract partition IDs
# K = length(part_id)                # number of partitions
# 
# 
# # ------------------------------------------------------------------------------
# 
# # exp_terms = vector("list", K)    # store terms in the exp in final approx
# candidate_psi = vector("list", K) # store the sub-partitions after resample
# optimal_psi = vector("list", K) # store the sub-partitions after resample
# 
# # k = 1
# 
# for (k in 1:K) {
#     
#     # print(k)
#     N_k_p = curr_part$n_obs[k] * n_samps  # num of samples to draw from A_k
#     
#     # set of lower/upper bounds corresponding to the k-th partition
#     # (1) omit the other three variables so that we only have ub/lb columns
#     # (2) transform row vector into a (D x 2) matrix of lb/ub
#     part_k = curr_part %>%                
#         dplyr::filter(leaf_id == part_id[k]) %>%
#         dplyr::select(-c(leaf_id, psi_star, n_obs)) %>% 
#         unlist() %>% matrix(ncol = 2, byrow = T)
#     
#     # sample unif from each lower/upper bound pair to form a D-dim vector
#     # lower bound in first column, upper bound in second column
#     # part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
#     
#     # sample uniformly from the D-dim partition whose lower/upper bounds
#     # are defined in part_k above
#     resamp_k = Matrix_runif(N_k_p, 
#                             lower = part_k[,1], 
#                             upper = part_k[,2]) %>% data.frame
#     
#     # compute psi(u) for each of the samples
#     u_df_k = preprocess(resamp_k, D, params) # N_k_p x (D_u + 1)
#     
#     ## NEW CODE 5/29 -------------------------------------------------------
#     
#     # optimal value of psi_star of the original k-th partition
#     c_k_star = curr_part$psi_star[k] 
#     
#     # compute R(u) for each of the samples; since we're in the k-th
#     # partition, we can directly subtract from psi_star[k], which we have
#     # computed and stored from the previous (original) call to hml_approx()
#     # this will need to be modified in the (to be)-implemented func that
#     # handles fitting the residuals, rather than the func values of psi(u)
#     R_df_k = u_df_k %>% dplyr::mutate(R_u = psi_u - c_k_star) %>% 
#         dplyr::select(-psi_u)
#     
#     # fit (u, R(u)) into decision tree to obtain partition
#     resid_tree = rpart(R_u ~ ., R_df_k)
#     
#     # extract the (data-defined) support from R_df_k
#     resid_support = extractSupport(R_df_k, D)
#     
#     # extract partition from resid_tree
#     # ntoe we omit 'psi_hat' the fitted value for the residual returned from 
#     # the tree, since this value results in an underestimation of the 
#     # log marginal likelihood. instead, use our own objective function to 
#     # choose the optimal value of R(u) for each (sub-)partition
#     resid_partition = extractPartition(resid_tree, resid_support) %>% 
#         dplyr::select(-psi_hat)
#     
#     ####  8/9 from here, we can just compute the log volume of each of the 
#     ####  hypercubes and store them separately (everything will be done within
#     ####  this function)
#     # part_volume = data.frame(leaf_id = resid_partition$leaf_id,
#     #                          log_vol = log_volume(resid_partition, D))
#     
#     # number of sub-partitions of partition k
#     s_k = nrow(resid_partition) 
#     
#     # compute opt value (chosen by tree) for each sub-partition
#     # e_kj_opt : leaf_id, Ru_choice, Ru_star, logJ_star, n_obs, lb/ub
#     
#     # 8/9 udpate: need two data structures
#     # (1) optimal psi_star that minimizes loss function
#     # (2) candidate psi_values to be used in exponential weighting
#     
#     e_kj = partition_opt_update(resid_tree, R_df_k, resid_partition, D)
#     # (1) obtain optimal: [leaf_id, psi_star]
#     opt = e_kj$Ru_df_opt
#     # (2) obtain candidates: [leaf_id, psi_1, ..., psi_5]
#     cand = e_kj$Ru_df_cand
#     
#     # compute psi_star = c_k_star + e_kj_star, e_kj_star = Ru_star
#     psi_star = cbind(leaf_id = opt$leaf_id, psi_star = opt$Ru_star + c_k_star,
#                      n_obs = opt$n_obs)
#     
#     # compute psi_tilde using the candidate R(u) values
#     psi_tilde = cbind(leaf_id = cand$leaf_id, cand[,-1] + c_k_star)
#     names(psi_tilde) = c("leaf_id", paste("psi_tilde_", 
#                                           c(1:(ncol(psi_tilde)-1)), sep = ''))
#     resid_partition = psi_star %>% merge(resid_partition, by = 'leaf_id')
#     
#     # append column of number of observations to resid_partitions
#     # sub_partitions[[k]] = resid_partition
#     
#     #### store (1) candidates, (2) optimal partition
#     candidate_psi[[k]] = psi_tilde
#     optimal_psi[[k]]   = resid_partition
#     
# } # end of for loop iterating over partitions
# 
# ## all candidates
# candidate_df = do.call(rbind, candidate_psi)
# 
# ## all optimal
# optimal_df = do.call(rbind, optimal_psi)
# 
# ## re-index the leaf IDs
# candidate_df = candidate_df %>% dplyr::mutate(leaf_id = 1:nrow(candidate_df))
# optimal_df = optimal_df %>% dplyr::mutate(leaf_id = 1:nrow(optimal_df))
# 
# candidate_df %>% head
# optimal_df %>% head
# 
# drops   = c('leaf_id', 'psi_star', 'n_obs')
# bounds  = optimal_df[,!(names(optimal_df) %in% drops)]
# log_vol = log_volume(optimal_df, drops, D)
# 
# exp_terms_mat = -candidate_df[,-1] + log_vol
# 
# apply(exp_terms_mat, 2, log_sum_exp)
# 
# #### determine which partitions the original samples lie in 
# u_sub = u_df %>% dplyr::select(-psi_u)
# 
# # extract only partition-defining columns
# part2 = bounds
# 
# # extract the optimal value for each partition
# psi_star_2 = candidate_df$psi_1
# 
# # convert partitions stored as K rows into a list of K partitions
# part2_list = lapply(split(part2, seq(NROW(part2))), matrix_part)
# 
# part2_id = apply(u_sub, 1, query_partition, part_list = part2_list)
# 
# u_df_wt = u_df_wt %>% dplyr::mutate(psi_2 = psi_star_2[part2_id])
# 
# # replace missing values with psi_1 values (this will be done for stage 3 too)
# u_df_wt = u_df_wt %>% dplyr::mutate(psi_2 = ifelse(is.na(psi_2), psi_1, psi_2))
# 
# u_df_wt$psi_2 %>% is.na %>% sum







































