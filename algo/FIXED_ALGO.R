



# ----------------------------------------------------------------------

#### at this point we're ready to compute the error associated with
#### original posterior samples using candidate psi values from BOTH
#### first and second stage partitioning

## compute approximations using first stage (original) samples

# optimal_df   : [leaf_id | psi_star | n_obs | u1_lb, u1_ub, ... uD_ub]
# candidate_df : [leaf_id | psi_2_1 ... psi_2_5 | psi_1_1 ... psi_1_5]
#
# psi_2_i, psi_1_i correspond to 2nd and 1st stage candidates, resp.
# optimal_df is used for its partition information to compute the vol
# of each partition
# candidate_df is used so that we know each partition's candidate psi 
# values
param = prior
set.seed(1)
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
psi_star = psi_update$psi_star
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

error = apply(psi_tilde_df, 2, logQ, u_df$psi_u) %>% 
    melt(value.name = 'error')

log_vol = log_volume(optimal_df, drops, D)
exp_terms_mat = -candidate_df[,-1] + log_vol
all_approx = reshape2::melt((apply(exp_terms_mat, 2, log_sum_exp)),
                              value.name = 'approx')

all_approx

approx_fs_error = merge(all_approx, error, by = 0)

approx_fs_error = approx_fs_error %>% dplyr::mutate(error = -log(J) + error)


approx_fs_error = approx_fs_error %>% 
    dplyr::mutate(wt_2d = exp(-error/(2*D))) %>% 
    dplyr::mutate(norm_wt_2d = wt_2d / sum(wt_2d)) %>% 
    dplyr::select(-c('wt_2d'))

approx_fs_error = approx_fs_error %>% 
    dplyr::mutate(wt_sqrt = exp(-error/(2*sqrt(D)))) %>% 
    dplyr::mutate(norm_wt_sqrt = wt_sqrt / sum(wt_sqrt)) %>% 
    dplyr::select(-c('wt_sqrt'))

approx_fs_error
(wt_approx1 = with(approx_fs_error, sum(approx * norm_wt_2d)))
(wt_approx2 = with(approx_fs_error, sum(approx * norm_wt_sqrt)))

# ------------------------------------------------------------------------------

num_ss = nrow(stage_s_part$u_df_resample)
psi_ss = stage_s_part$u_df_resample[,-c(1:D)]
# check the 1st stage candidates have the same number of unique values
table(psi_ss$psi_2_2) %>% length 
# check the number of 1st stage partitions == number of 1st stage values
hml_approx$n_partitions

# compute ss error
ss_error = apply(psi_ss %>% dplyr::select(-psi_u), 2, logQ, psi_ss$psi_u) %>% 
    melt(value.name = 'ss_error')
total_error = error %>% merge(ss_error, by = 0)
total_error = total_error %>% dplyr::mutate(ss_error = - log(num_ss) + ss_error)
names(total_error)= c("psi_candidates", "fs_error", "ss_error")

total_error = total_error %>% mutate(total_error = fs_error + ss_error)


approx_tot_error =  merge(all_approx %>% 
                              dplyr::mutate(psi_candidates = row.names(all_approx)), 
                          total_error, by = 'psi_candidates')


approx_tot_error = approx_tot_error %>% 
    dplyr::mutate(wt_2d = exp(-total_error/(2*D))) %>% 
    dplyr::mutate(norm_wt_2d = wt_2d / sum(wt_2d)) %>% 
    dplyr::select(-c('wt_2d'))

approx_tot_error = approx_tot_error %>% 
    dplyr::mutate(wt_sqrt = exp(-total_error/(2*sqrt(D)))) %>% 
    dplyr::mutate(norm_wt_sqrt = wt_sqrt / sum(wt_sqrt)) %>% 
    dplyr::select(-c('wt_sqrt'))

approx_tot_error
(wt_approx1 = with(approx_tot_error, sum(approx * norm_wt_2d)))
(wt_approx2 = with(approx_tot_error, sum(approx * norm_wt_sqrt)))


# ------------------------------------------------------------------------------

all_df = rbind(psi_tilde_df, psi_ss %>% dplyr::select(-psi_u))
all_psi = c(u_df$psi_u, psi_ss$psi_u)

all_error = (-log(J + num_ss) + apply(all_df, 2, logQ, all_psi)) %>% 
    melt(value.name = 'error')

approx_all_error =  merge(all_approx, all_error, by = 0)

approx_all_error = approx_all_error %>% 
    dplyr::mutate(wt_2d = exp(-error/(2*D))) %>% 
    dplyr::mutate(norm_wt_2d = wt_2d / sum(wt_2d)) %>% 
    dplyr::select(-c('wt_2d'))

approx_all_error = approx_all_error %>% 
    dplyr::mutate(wt_sqrt = exp(-error/(2*sqrt(D)))) %>% 
    dplyr::mutate(norm_wt_sqrt = wt_sqrt / sum(wt_sqrt)) %>% 
    dplyr::select(-c('wt_sqrt'))

approx_all_error
(wt_approx1 = with(approx_all_error, sum(approx * norm_wt_2d)))
(wt_approx2 = with(approx_all_error, sum(approx * norm_wt_sqrt)))

