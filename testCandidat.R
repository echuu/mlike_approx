



## iteratively make columns to populate the psi values for each row

# extract posterior samples without Psi(u) evaluation
u_sub = u_df %>% dplyr::select(-psi_u)
#### PREPARE FOR STAGE 1/2/3 PROCESSING ----------------------------------------


# u_df_star: [psi_1_star, psi_2_star, psi_3_star | u, psi_u]
# contains each stage's optimal value of psi for each point 
# this lets us directly access the *previous* stage's psi_star values

# u_df_star = data.frame(matrix(0, J, 3), u_df) # don't need u_df part of it
psi_star = data.frame(matrix(0, J, 3)) # don't need u_df part of it
psi_star %>% head

all_psi_candidates = data.frame(matrix(0, J, 3 * 5)) # num_stages * num_cands

# extract only partition-defining columns
#### STAGE 1 -------------------------------------------------------------------

stage = 1
part1 = og_part %>% dplyr::select(-c('leaf_id', 'psi_star', 'n_obs'))

# extract the optimal value for each partition
psi_star_1 = og_part$psi_star

# convert partitions stored as K rows into a list of K partitions
part1_list = lapply(split(part1,seq(NROW(part1))), matrix_part)
part1_id = apply(u_sub, 1, query_partition, part_list = part1_list)
part1_id %>% head

u_df_wt = u_df %>% dplyr::mutate(psi_1 = psi_star_1[part1_id])

psi_star[,stage] = psi_star_1[part1_id]

psi_star %>% head


#### STAGE 2 -------------------------------------------------------------------

stage1_part = og_part
stage2_fit = next_stage(stage1_part, D, n_samps, prior)

part2_cand = stage2_fit$candidate_psi
part2_opt  = stage2_fit$optimal_part

part2_cand %>% head
part2_opt %>% head


psi2_update = compute_weights(u_df, stage2_fit, 2, psi_star)
psi3_update$approx_error
psi_star = psi3_update$psi_star


#### STAGE 3 -------------------------------------------------------------------


stage2_part = stage2_fit$optimal_part
stage3_fit = next_stage(stage2_part, D, n_samps / 2, prior)

part3_cand = stage3_fit$candidate_psi
part3_opt  = stage3_fit$optimal_part

part3_cand %>% head
part3_opt %>% head


psi_update = compute_weights(u_df, stage3_fit, 3, psi_star)
psi_update$approx_error
psi_update$psi_star %>% head



#### finish running hml_const() function

params = prior
approx = vector("list", 3)

prev_part = og_part
n_samps = 10
for (s in 2:3) {
    
    stage_s_part = next_stage(prev_part, D, n_samps, params)
    psi_update = compute_weights(u_df, stage_s_part, 2, psi_star)
    
    approx[[s]] = psi_update$approx_error
    
    prev_part = stage_s_part$optimal_part
    psi_star = psi_update$psi_star # used in next iteration
    n_samps = n_samps / 2
}

approx[[3]]

do.call(rbind, approx)




all_psi_candidates = data.frame(matrix(0, J, 3 * 5)) # num_stages * num_cands







part2 = part2_opt[,!(names(part2_opt) %in% drops)]

# convert partitions stored as K rows into a list of K partitions
part2_list = lapply(split(part2, seq(NROW(part2))), matrix_part)

## these IDs are the SAME for all candidates, so we use these to index through 
## all the candidate values
part2_id = apply(u_sub, 1, query_partition, part_list = part2_list)

# points that fall outside re-sampled partitions - in this case, use previous
# stage's points 
use_prev_id = which(is.na(part2_id))

## go through each candidate
# create new dataframe to store the values (rows match u_df_wt)
n_cand = ncol(candidate_df) - 1 # subtract off the leaf_id column
psi_tilde_df = data.frame(matrix(0, J, n_cand))
stage = 2
for (i in 1:n_cand) {
    
    # identify the psi candidate value from candidate_df
    psi_cand_i = candidate_df[,i+1]
    psi_tilde_df[,i] = psi_cand_i[part2_id]
    psi_tilde_df[use_prev_id,i] = psi_star[use_prev_id, stage-1]
}

psi_star[,stage] = part2_opt$psi_star[part2_id]
psi_star[,stage][use_prev_id] = psi_star[use_prev_id, stage-1]
psi_star %>% head

psi_tilde_df %>% head




## populate the 'all_psi_candidates' dataframe with these candidate values
start = (stage-1)*n_cand+1
end = start + n_cand - 1
all_psi_candidates[,start:end] = psi_tilde_df


#### IN PROGRESS ---------------------------------------------------------------

#### STAGE 3 -------------------------------------------------------------------


stage2_part = stage2_fit$optimal_part
stage3_fit = next_stage(stage2_part, D, n_samps / 2, prior)

part3_cand = stage3_fit$candidate_psi
part3_opt  = stage3_fit$optimal_part

part3_cand %>% dim
part3_opt %>% dim


psi_update = compute_weights(u_df, stage3_fit, 3, psi_star)
psi_update$approx_error
psi_update$psi_star %>% head

compute_approx(stage3_fit)


######### compute the log ML approx using all 5 candidates



######### this part is only used for computing the errors, weights

part3 = optimal_df[,!(names(optimal_df) %in% drops)]
part3 %>% head
# convert partitions stored as K rows into a list of K partitions
part3_list = lapply(split(part3, seq(NROW(part3))), matrix_part)

## these IDs are the SAME for all candidates, so we use these to index through 
## all the candidate values
part3_id = apply(u_sub, 1, query_partition, part_list = part3_list)

# points that fall outside re-sampled partitions - in this case, use previous
# stage's points 
use_prev_id = which(is.na(part3_id))
use_prev_id %>% length

## go through each candidate
# create new dataframe to store the values (rows match u_df_wt)
n_cand = ncol(candidate_df) - 1 # subtract off the leaf_id column
psi_tilde_df = data.frame(matrix(0, J, n_cand))
stage = 3
for (i in 1:n_cand) {
    
    # identify the psi candidate value from candidate_df
    psi_cand_i = candidate_df[,i+1]
    psi_tilde_df[,i] = psi_cand_i[part2_id]
    psi_tilde_df[use_prev_id,i] = psi_star[use_prev_id, stage-1]
}

psi_star[,stage] = part3_opt$psi_star[part3_id]
psi_star[,stage] %>% is.na %>% sum
psi_star[,stage][use_prev_id] = psi_star[use_prev_id, stage-1]
psi_star[,stage] %>% is.na %>% sum
psi_star %>% head


## populate the 'all_psi_candidates' dataframe with these candidate values
start = (stage-1)*n_cand+1
end = start + n_cand - 1
all_psi_candidates[,start:end] = psi_tilde_df

psi_tilde_df %>% head
psi_tilde_df %>% dim



drops   = c('leaf_id', 'psi_star', 'n_obs')
bounds  = optimal_df[,!(names(optimal_df) %in% drops)]
log_vol = log_volume(optimal_df, drops, D)

exp_terms_mat = -candidate_df[,-1] + log_vol

exp_terms_mat %>% dim

apply(exp_terms_mat, 2, log_sum_exp)


# obtain error for each of the candidates
apply(all_psi_candidates, 2, MAE, u_df_wt$psi_u)


all_psi_candidates %>% head

#### append all candidates to a master candidate df




