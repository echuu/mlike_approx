



c_k = max(u_df$psi_u)
c_star = min(u_df$psi_u)

c_k = 2
c_star = 3

abs(exp(-c_k) - exp(-c_star)) / exp(-c_k)

1 - exp(-c_star + c_k)

# this can be passed into log-sum-exp as one of the elements of the summation
# in the case that c_star > c_k
log1mexp(c_star - c_k) 

# in the case that c_k > c_star
log(abs(1 - exp(-c_star + c_k)))

log(abs((exp(-c_k) - exp(-c_star))) / exp(-c_k))

(log_kj = log1mexp(c_k - c_star) - c_star + c_k) # answer

log_sum_exp(log_kj)

log(abs(1 - exp(-c_star + c_k)))



# compute the log-loss function for a given c_star and partition A_k

# leaf_id, psi_hat
psi_hat_df = hml_approx$param_out %>% dplyr::select(leaf_id, psi_hat)

u_df_all = hml_approx$u_df_fit
u_df_all = u_df_all %>% dplyr::select(leaf_id, psi_u, psi_star)
u_df_all %>% head

psi_all = u_df_all %>% dplyr::select(leaf_id, psi_u, psi_star)
psi_all = merge(psi_all, psi_hat_df, by = 'leaf_id')

# psi_all: leaf_id, psi_u, psi_star, psi_hat

part_k = psi_all %>% dplyr::filter(leaf_id == 14)

# compute log(Q(c_star)) for part_k


c_k = part_k$psi_u               # true values of psi for partition_k
c_star = part_k$psi_hat[1]       # choice 1: fitted value from tree
c_star = part_k$psi_star[1]      # choice 2: max value of psi over partition

c_star = part_k$psi_u %>% max    # choice 3: median value of psi over partition
c_star = part_k$psi_u %>% median # choice 3: median value of psi over partition
c_star = part_k$psi_u %>% mean   # choice 4: mean value of psi over partition


L_k = nrow(part_k) # can use n_obs later
log_rel_error = rep(NA, L_k) # store the log relative error

# l_k = 6
for (l_k in 1:L_k) {
    ## perform a stable calculation of log(abs(1-exp(-c_star + c_k[l_k])))
    ## by considering cases when c_star > c_k[l_k], c_star < c_k[l_k] 
    
    if (c_star > c_k[l_k])
        log_rel_error[l_k] = log1mexp(c_star - c_k[l_k])
    else if (c_star < c_k[l_k])
        log_rel_error[l_k] = log1mexp(c_k[l_k] - c_star) - c_star + c_k[l_k]
    
    # if c_star == c_k[l_k] : do nothing; NA value will be skipped over in 
    # final calculation
    
} # end of for loop iterating over each element in k-th partition


# log_rel_error_og = rep(NA, L_k)
# for (l_k in 1:L_k) {
#     ## original calculation w/o stabilizing the computation
#     if (c_star != c_k[l_k])
#         log_rel_error_og[l_k] = log(abs(1-exp(-c_star + c_k[l_k])))
# } # end of for loop iterating over each element in k-th partition


(logQ_psi_hat = log_sum_exp(log_rel_error[!is.na(log_rel_error)])) # 5.422897


psi_all = hml_approx$u_df_fit %>% dplyr::select(leaf_id, psi_u, psi_star)
c_k = head(psi_all$psi_u, 100)    # true values of psi for partition_k
c_star = min(c_k)                 # choice 2: max value of psi over partition

logQ(c_star, c_k)
logQ_fast(c_star, c_k)

x = runif(100)
stopifnot(all.equal(logQ(c_star, x), logQ_fast(c_star, x)))

library(microbenchmark)
microbenchmark(
    logQ(c_star, c_k),
    logQ_fast(c_star, c_k)
)


x = runif(100)
stopifnot(all.equal(logJ(c_star, x), logJ_fast(c_star, x)))

microbenchmark(
    logJ(c_star, c_k),
    logJ_fast(c_star, c_k)
)


logJ(c_star, c_k)
logJ_fast(c_star, c_k)


logJ_fast = function(c_star, c_k) {
    
    L = length(c_k)
    
    log_rel_error = rep(NA, L) # store the log relative error
    
    # for (l in 1:L) {
    #     ## perform a stable calculation of log(abs(1-exp(-c_star + c_k[l])))
    #     ## by considering cases when c_star > c_k[l], c_star < c_k[l] 
    #     if (c_star < c_k[l])
    #         log_rel_error[l] = log1mexp(c_k[l] - c_star)
    #     else if (c_star > c_k[l])
    #         log_rel_error[l] = log1mexp(c_star - c_k[l]) + c_star - c_k[l]
    #     
    #     # if c_star == c_k[l_k] : do nothing; NA value will be skipped over in 
    #     # final calculation
    #     
    # } # end of for loop iterating over each element in k-th partition
    
    sign_ind = (c_star > c_k)
    log_rel_error = log1mexp(abs(c_star - c_k)) + (c_star - c_k) * sign_ind
    
    return(log_sum_exp(log_rel_error[!is.na(log_rel_error)]))
    
}





## logQ() function -------------------------------------------------------------
# c_star : value whose cost function is to be evaluated
# c_k    : L-dim vector of function evaluations in the k-th partition
logQ_fast = function(c_star, c_k) {
    
    L = length(c_k)
    
    log_rel_error = rep(NA, L) # store the log relative error
    
    # for (l in 1:L) {
    #     ## perform a stable calculation of log(abs(1-exp(-c_star + c_k[l])))
    #     ## by considering cases when c_star > c_k[l], c_star < c_k[l] 
    #     if (c_star > c_k[l]) {
    #         log_rel_error[l] = log1mexp(c_star - c_k[l])
    #     } else if (c_star < c_k[l]) {
    #         log_rel_error[l] = log1mexp(c_k[l] - c_star) - c_star + c_k[l]
    #     }
    #     
    #     # if c_star == c_k[l_k] : do nothing; NA value will be skipped over in 
    #     # final calculation
    #     
    # } # end of for loop iterating over each element in k-th partition
    
    sign_ind = (c_star < c_k)
    log_rel_error = log1mexp(abs(c_star - c_k)) + (c_k - c_star) * sign_ind
    
    return(log_sum_exp(log_rel_error[!is.infinite(log_rel_error)]))
    
}
# end of logQ() function -------------------------------------------------------




## logMSE() function -------------------------------------------------------------
# c_star : value whose cost function is to be evaluated
# c_k    : L-dim vector of function evaluations in the k-th partition
logMSE = function(c_star, c_k) {
    
    L = length(c_k)
    
    log_mse = rep(NA, L) # store the log mse
    
    for (l in 1:L) {
        if (c_star > c_k[l]) {
            log_mse[l] = 2 * (-c_k[l] + log1mexp(c_star - c_k[l]))
        } else if (c_star < c_k[l]) {
            log_mse[l] = 2 * (-c_star + log1mexp(c_k[l] - c_star))
        }
    }
        
    return(- log(L) + log_sum_exp(log_mse[!is.na(log_mse)]))
    
}
# end logMSE() function --------------------------------------------------------
















