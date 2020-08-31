



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

# sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
# beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
# u_samp = data.frame(beta_mat, sigmasq_post)
# u_df = preprocess(u_samp, D, prior)

source("setup.R")
param = prior
# param = params

set.seed(1)
# D = D_u
n_stage = 2
n_samps = 10
params = param
approx = vector("list", n_stage)
psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it

logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
prev_part = logml_approx$param_out$optimal_part # first stage partition
# psi_update = compute_weights(u_df, logml_approx$param_out, 1, psi_star)

# --------------------------------------------------------------------------

n_cand = 8
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
stage_s_part$u_df_resample %>% dim
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



# all_error = (-log(J + num_ss) + apply(all_df, 2, logQ, all_psi)) %>%
#     melt(value.name = 'error')
all_error = (-log(J + num_ss) + (apply(all_df, 2, logJ, all_psi) + 
                                     apply(all_df, 2, logQ, all_psi))/2) %>% 
    melt(value.name = 'error')

# apply(all_df[,1:5],  2, logJ, all_psi)
# apply(all_df[,6:10], 2, logQ, all_psi)


all_error = rbind(apply(all_df[,1:8], 2,  logQ, all_psi) %>%
                      melt(value.name = 'error'),
                  apply(all_df[,9:16], 2, logJ, all_psi) %>%
                      melt(value.name = 'error')) -log(J + num_ss)

approx_all_error =  merge(all_approx, all_error, by = 0)

approx_all_error = approx_all_error %>% 
    dplyr::mutate(wt_1 = exp(-error/(6*D))) %>% 
    dplyr::mutate(norm_wt_1 = wt_1 / sum(wt_1)) %>% 
    dplyr::select(-c('wt_1'))

approx_all_error = approx_all_error %>% 
    dplyr::mutate(wt_2 = exp(-error/(8*D))) %>% 
    dplyr::mutate(norm_wt_2 = wt_2 / sum(wt_2)) %>% 
    dplyr::select(-c('wt_2'))

approx_all_error = approx_all_error %>% 
    dplyr::mutate(wt_3 = exp(-error/(100*D))) %>% 
    dplyr::mutate(norm_wt_3 = wt_3 / sum(wt_3)) %>% 
    dplyr::select(-c('wt_3'))

approx_all_error
(wt_approx1 = with(approx_all_error, sum(approx * norm_wt_1)))
(wt_approx2 = with(approx_all_error, sum(approx * norm_wt_2)))
(wt_approx3 = with(approx_all_error, sum(approx * norm_wt_3)))







