
param = prior

curr_part = prev_part
drops = c('leaf_id', 'psi_star', 'n_obs')


#### compute all of the resampled error ----------------------------------------

# this has stage 1 error for 2nd stage samples
resamp_error = function(part_fit) {
    psi_tilde_df = part_fit$u_df_resample %>% 
        dplyr::select(psi_tilde_1:psi_1)
    error = apply(psi_tilde_df, 2, logQ, part_fit$u_df_resample$psi_u)
    return(error)
}

psi_update = compute_weights(u_df, stage_s_part, s, psi_star)
psi_update$approx_error
psi_update$psi_tilde_df %>% head

ss_error = melt(resamp_error(stage_s_part), value.name = 'error')

approx_error = psi_update$approx_error %>% 
    dplyr::mutate(ss_error = ss_error$error) %>% 
    dplyr::mutate(tot_error = log(exp(error) + exp(ss_error)))

# create a version of u_df where each row is a re-sampled point with columns:
# (0) actual psi values for each of the points
# (1) psi values from first stage candidates  (5 candidates)
# (2) psi values from second stage candidates (5 candidates)

# input
stage_s_part$u_df_resample %>% head

stage_s_part$u_df_resample %>% head



curr_part = prev_part
psi_candidates
n_cand = 5
params = prior

part_id = curr_part$leaf_id        # extract partition IDs
K = length(part_id)                # number of partitions

candidate_psi = vector("list", K) # store the sub-partitions after resample
opt_psi_part  = vector("list", K) # store the sub-partitions after resample
u_df_k_list   = vector("list", K) # store resampled points w/ psi, psi_star

N_k_p = curr_part$n_obs[k] * n_samps  # num of samples to draw from A_k

# set of lower/upper bounds corresponding to the k-th partition
# (i)  omit the other three variables so that we only have ub/lb columns
# (ii) transform row vector into a (D x 2) matrix of lb/ub
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

#### (1) psi values from first stage candidates  (5 candidates)
# save candidate values for this partition by extracting it from psi_candidates
psi_cand_k = psi_candidates %>% dplyr::filter(leaf_id == part_id[k]) %>% 
    dplyr::select(-c('leaf_id'))
names(psi_cand_k) = paste("psi_1", 1:n_cand, sep = '_')
u_df_k = u_df_k %>% dplyr::mutate(psi_cand_k)

u_df_k_cand = cbind(u_df_k, psi_cand_k)


#### (2) psi values from second stage candidates (5 candidates)






















#### TODO: compute first stage errors  -----------------------------------------

u_df_resample = stage_s_part$u_df_resample
psi_tilde_df = u_df_resample %>% dplyr::select(psi_tilde_1:psi_tilde_5)

error_resample = apply(psi_tilde_df, 2, logQ, u_df_resample$psi_u)

compute_weights(u_df, stage_s_part, s, psi_star)

part_fit = stage_s_part
s = 2

u_sub = u_df %>% dplyr::select(-psi_u)


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


