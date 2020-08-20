

## compute_expterms() ----------------------------------------------------------
# part_info: df containing leaf_id, psi_star, lb/ub of D-intervals that make 
#            up each of the K partitions
# D        : dimension of parameter
log_part_volume = function(bounds, D) {
    
    log_vol = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
        log() %>% rowSums()
    
    # (1) compute length of each of the D intervals for each partition
    # (2) take log of each interval difference
    # (3) add up each row for log volume of each partition
    log_vol = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>% 
        log() %>% rowSums()
    
    return(log_vol) # K-dim vector
} # end of compute_expterms() function -----------------------------------------


n_stage = 2
n_samps = 10
params = prior
param = prior

# set.seed(1)
# approx = vector("list", n_stage)
# psi_star = data.frame(matrix(0, J, n_stage)) # don't need u_df part of it

logml_approx = hybrid_ml(D, u_df, J, param)     # fit first stage

prev_part = logml_approx$param_out$optimal_part # first stage partition
psi_candidates = logml_approx$param_out$candidate_psi
names(psi_candidates) = c('leaf_id', paste("psi_1", 1:n_cand, sep = '_'))

## create psi_star column
u_df_star = u_df %>% 
    dplyr::mutate(leaf_id = logml_approx$u_rpart$where) %>% 
    plyr::join(prev_part %>% dplyr::select(c('leaf_id', 'psi_star')),
               by = 'leaf_id') %>% 
    dplyr::select(-c('leaf_id'))


##### perform second stage -----------------------------------------------------
curr_part = prev_part
n_cand = 5


part_id = curr_part$leaf_id        # extract partition IDs
K = length(part_id)                # number of partitions

candidate_psi = vector("list", K) # store the sub-partitions after resample
opt_psi_part  = vector("list", K) # store the sub-partitions after resample
u_df_k_list   = vector("list", K) # store resampled points w/ psi, psi_star

for (k in 1:K) {
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
    # save candidate values for this partition by extracting from psi_candidates
    psi_cand_k = psi_candidates %>% dplyr::filter(leaf_id == part_id[k]) %>% 
        dplyr::select(-c('leaf_id'))
    # names(psi_cand_k) = paste("psi_1", 1:n_cand, sep = '_')
    # u_df_k = u_df_k %>% dplyr::mutate(psi_cand_k)
    
    # u_df_k_cand = data.frame(u_df_k, psi_cand_k)
    # u_df_k_cand %>% head
    
    # verify order is preserved
    # (u_df_k_cand$psi_u != u_df_k$psi_u) %>% sum
    
    #### (2) psi values from second stage candidates (5 candidates)
    
    #### (2a) do we even need the optimal value if we're doing candidates? I 
    ####  don't think so -- just create another candidate matrix. At this point, 
    ####  the only optimal value that we rely on the c_k_star from the PREVIOUS 
    ####  (first) stage so that we can form the candidates, i.e., 
    ####  c_k_star + (cand1, cand2, ..., cand_5)
    
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
    
    e_kj = partition_opt_update(resid_tree, R_df_k, resid_partition, D)
    # (1) obtain optimal: [leaf_id, psi_star]
    opt = e_kj$Ru_df_opt
    # (2) obtain candidates: [leaf_id, psi_1, ..., psi_5]
    cand = e_kj$Ru_df_cand
    
    #### the following calculation is superfluous since we don't go beyond 
    #### 2nd stage but keep it for now since it's low-cost and easy to compute 
    #### (might need it later at some point)
    # compute psi_star = c_k_star + e_kj_star, e_kj_star = Ru_star
    psi_star = cbind(leaf_id  = opt$leaf_id, 
                     psi_star = opt$Ru_star + c_k_star,
                     n_obs    = opt$n_obs)
    
    resid_partition = psi_star %>% merge(resid_partition, by = 'leaf_id')
    
    
    # store u_df_k with psi_star value so that we can compute the errors
    # u_df_k = u_df_k %>% dplyr::mutate(leaf_id = resid_tree$where)
    
    u_df_k_cand = u_df_k %>% dplyr::mutate(leaf_id = resid_tree$where)
    
    # compute psi_tilde using the candidate R(u) values
    # first column is leaf_id, throw this out so we can do element-
    # wise multiplication
    # psi_tilde is used later to compute the final approximation
    psi_tilde = cbind(leaf_id = cand$leaf_id, cand[,-1] + c_k_star)
    names(psi_tilde) = c("leaf_id", 
                         paste("psi_2_", c(1:(ncol(psi_tilde)-1)), sep = ''))
    # why do we combine this with the first stage candidates??
    # merge so that u_df_k_cand has access to previous stage psi candidates
    psi_tilde = data.frame(psi_tilde, psi_cand_k)
    
    u_df_k_cand = merge(u_df_k_cand, psi_tilde, by = 'leaf_id')
    u_df_k_cand %>% head
    
    candidate_psi[[k]] = psi_tilde
    opt_psi_part[[k]]  = resid_partition
    u_df_k_list[[k]]   = u_df_k_cand %>% dplyr::select(-c('leaf_id')) 
}

### psi_tilde       -> used to compute final approximation using each candidate
### u_df_k_cand     -> used to compute the error + weights for each logml approx
### resid_partition -> used to compute the volume of each of the partitions

## all candidates
candidate_df = do.call(rbind, candidate_psi)

## all optimal
optimal_df = do.call(rbind, opt_psi_part)

## all of u_df_k 
u_df_resample = do.call(rbind, u_df_k_list)


## re-index the leaf IDs
K = nrow(candidate_df)
candidate_df = candidate_df %>% dplyr::mutate(leaf_id = 1:K)
optimal_df = optimal_df %>% dplyr::mutate(leaf_id = 1:K)
##### end second stage ---------------------------------------------------------

# ss_out = next_stage(prev_part, psi_candidates, D, n_samps, params)
# andidate_df = ss_out$candidate_psi
# names(candidate_df) = c("leaf_id", 
#                         paste("psi_2_", c(1:(ncol(candidate_df)-1)), sep = ''))
# optimal_df = ss_out$optimal_part
# u_df_resample = ss_out$u_df_resample


### psi_tilde       -> used to compute final approximation using each candidate
### u_df_k_cand     -> used to compute the error + weights for each logml approx
### resid_partition -> used to compute the volume of each of the partitions

# 2nd stage candidates
candidate_df %>% head
candidate_df %>% dim

optimal_df %>% head
u_df_resample %>% head


###### compute errors all together after both stages are finished

#### compute first stage errors -------------------------------------------------

#### add first stage candidate columns
## create candidate columns
u_df_cand = u_df %>% 
    dplyr::mutate(leaf_id = logml_approx$u_rpart$where) %>% 
    plyr::join(psi_candidates, by = 'leaf_id') %>% 
    dplyr::select(-c('leaf_id'))

#### add second stage candidate columns
u_sub = u_df %>% dplyr::select(-psi_u)
drops   = c('leaf_id', 'psi_star', 'n_obs')

# partitions from 2nd stage partition
part = optimal_df[,!(names(optimal_df) %in% drops)]
part_list = lapply(split(part, seq(NROW(part))), matrix_part)

# for each posterior sample, determine which partition it belongs in
part_id = apply(u_sub, 1, query_partition, part_list = part_list)

# for posterior samples that fall outside of the re-sampled partitions,
# we will use the psi_star value obtained during the previous stage
use_prev_id = which(is.na(part_id))
stage_1_psi_star = u_df_star$psi_star[use_prev_id]

# n_cand = ncol(candidate_df) - 1 # subtract off the leaf_id column
n_cand = 5 # subtract off the leaf_id column
psi_tilde_df = data.frame(matrix(0, nrow(u_df), n_cand))
names(psi_tilde_df) = names(candidate_df)[2:(n_cand+1)]
for (i in 1:n_cand) {
    # identify the psi candidate value from candidate_df
    psi_cand_i = candidate_df[,i+1]
    # match these with the corresponding partition ID for each point
    psi_tilde_df[,i] = psi_cand_i[part_id]
    # for points that fall outside of a 2nd stage part, use 1st stage part
    psi_tilde_df[use_prev_id,i] = stage_1_psi_star
}

psi_tilde_df %>% head
u_df_cand %>% head

all_psi_fs = cbind(psi_tilde_df, 
                   u_df_cand %>% dplyr::select(-c(paste('u', 1:D, sep = ''))))

fs_error = apply(all_psi_fs %>% dplyr::select(-c('psi_u')), 
                 2, logQ, all_psi_fs$psi_u)

#### compute second stage errors -----------------------------------------------

all_psi_ss = u_df_resample %>% dplyr::select(-c(paste('u', 1:D, sep = '')))
ss_error = apply(all_psi_ss %>% dplyr::select(-c(psi_u)), 
                 2, logQ, all_psi_ss$psi_u)

total_error = log(exp(melt(fs_error, value.name = 'error')) + 
                  exp(melt(ss_error, value.name = 'error')))


#### compute approximation for each of the candidates

## compute first stage approximation (5 candidates)
fs_part = prev_part %>% dplyr::select('u1_lb':paste('u', D, '_ub', sep = ''))
log_vol_1 = log_part_volume(fs_part, D)
exp_terms_mat_1 = -(psi_candidates %>% dplyr::select(-c('leaf_id'))) + log_vol_1
logml_approx_1 = reshape2::melt((apply(exp_terms_mat_1, 2, log_sum_exp)),
                              value.name = 'approx')


## compute second stage approximation
candidate_df %>% dim
part %>% dim

ss_part = part
log_vol_2 = log_part_volume(part, D)
exp_terms_mat_2 = -candidate_df[,2:(n_cand+1)] + log_vol_2
# test = -ss_out$candidate_psi[,-1] + log_vol_2
logml_approx_2 = reshape2::melt((apply(exp_terms_mat_2, 2, log_sum_exp)),
                              value.name = 'approx')

cbind(rbind(logml_approx_2, logml_approx_1), total_error)
      
all_approx = rbind(logml_approx_2, logml_approx_1)
all_approx = merge(all_approx, total_error, by = 0)

all_approx = all_approx %>% dplyr::mutate(wt_2d = exp(-error/(2*D))) %>% 
    dplyr::mutate(norm_wt_2d = wt_2d / sum(wt_2d)) %>% 
    dplyr::select(-c('wt_2d'))

all_approx = all_approx %>% dplyr::mutate(wt_sqrt = exp(-error/(2*sqrt(D)))) %>% 
    dplyr::mutate(norm_wt_sqrt = wt_sqrt / sum(wt_sqrt)) %>% 
    dplyr::select(-c('wt_sqrt'))

all_approx
(wt_approx1 = with(all_approx, sum(approx * norm_wt_2d)))
(wt_approx2 = with(all_approx, sum(approx * norm_wt_sqrt)))




