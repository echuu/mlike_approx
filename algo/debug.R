



param = prior
n_cand = 5

set.seed(1)
n_stage = 2
n_samps = 10
params = param
approx = vector("list", n_stage)
psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it

logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage
prev_part = logml_approx$param_out$optimal_part # first stage partition

# new code ---------------------------------------------------------------------
psi_candidates = logml_approx$param_out$candidate_psi
names(psi_candidates) = c('leaf_id', paste("psi_1", 1:n_cand, sep = '_'))

## create psi_star column
u_df_star = u_df %>% 
    dplyr::mutate(leaf_id = logml_approx$u_rpart$where) %>% 
    plyr::join(prev_part %>% dplyr::select(c('leaf_id', 'psi_star')),
               by = 'leaf_id') %>% 
    dplyr::select(-c('leaf_id'))


# ------------------------------------------------------------------------------


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

all_approx
wt_approx1
wt_approx2
