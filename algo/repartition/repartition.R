

library(rBayesianOptimization)


# contains: leaf_id, u1_lb, u1_ub, ... , uD_lb, uD_ub, n_obs
part_0 = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

K = length(part_set)

k = 1




N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from partition k
part_k = part_0 %>%           # set of lower/upper bounds
    dplyr::filter(leaf_id == part_set[k]) %>% 
    dplyr::select(-c(leaf_id, n_obs))

# sample uniformly from each lower/upper bound pair to form a D-dim vector
part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
apply(part_k_long, 1, runif)

runif(1, min = part_k_long[,1], max = part_k_long[,2])

resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1], 
                        upper = part_k_long[,2]) %>% data.frame

u_df_k = preprocess(resamp_k, D_u, param_list) # N_k_p x (D_u + 1)

c_k_approx = hml_const(1, D_u, u_df, N_k_p, param_list)

c_k_approx$const_approx # these are the elements in the exponential that we care about


# ------------------------------------------------------------------------------

# contains: leaf_id, u1_lb, u1_ub, ... , uD_lb, uD_ub, n_obs
part_0 = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

orig_partition = hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs))


K = length(part_set)

# initialize a list to store the vector containing the terms in exponential
# for each of the sub-partitions
# kth elmt is an s_k dim vector of terms that are to be exponentiated
# at the very end, all entries are unlisted and evaluated via log-sum-exp
exp_terms = vector("list", K) 
ck_star_list = vector("list", K)

for (k in 1:K) {
    
    
    PERC_K = orig_partition[k,]$perc
    
    if (PERC_K > 0.2) {
        print("using original partition")
        exp_terms[[k]] = hml_approx$const_approx[k]
        
        # N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        # part_k = part_0 %>%           # set of lower/upper bounds
        #     dplyr::filter(leaf_id == part_set[k]) %>%
        #     dplyr::select(-c(leaf_id, n_obs))
        # 
        # # sample uniformly from each lower/upper bound pair to form a D-dim vector
        # part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        # 
        # resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1],
        #                         upper = part_k_long[,2]) %>% data.frame
        # 
        # u_df_k = preprocess(resamp_k, D_u, param_list) # N_k_p x (D_u + 1)
        # 
        # c_k_approx = hml_const_mod(1, D_u, u_df_k, N_k_p, param_list)
        # 
        # ck_star_list[[k]] = c_k_approx$param_out %>%
        #     dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)
        # 
        # exp_terms[[k]] = c_k_approx$const_approx
        
    } else {
        N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        part_k = part_0 %>%           # set of lower/upper bounds
            dplyr::filter(leaf_id == part_set[k]) %>% 
            dplyr::select(-c(leaf_id, n_obs))
        
        # sample uniformly from each lower/upper bound pair to form a D-dim vector
        part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        
        resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1], 
                                upper = part_k_long[,2]) %>% data.frame
        
        u_df_k = preprocess(resamp_k, D_u, param_list) # N_k_p x (D_u + 1)
        
        c_k_approx = hml_const(1, D_u, u_df_k, N_k_p, param_list)
        
        ck_star_list[[k]] = c_k_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) 
        
        exp_terms[[k]] = c_k_approx$const_approx
    }
    
}


exp_terms
all_terms = exp_terms %>% unlist
all_terms


# -1183.447
# -1319.849
log_sum_exp(all_terms) # using "second-stage partitioning"
(approx_logml = hml_approx$const_vec) # using original partitioning scheme
(true_logml = lil(param_list))        # true log ML 


rpart_k = rpart(psi_u ~., u_df_k, cp = 0.0001)
rpart_k = rpart(psi_u ~., u_df_k, cp = 0.05)
plot(rpart_k)





#### hml_const() ---------------------------------------------------------------
#
#
hml_const_mod = function(N_approx, D, u_df_full, J, prior) {
    
    const_vec  = numeric(N_approx) # store constant approximation
    
    # compute approximation to LIL N_approx times
    for (t in 1:N_approx) {
        
        ## (1) subset out rows in u_df_full to be used in the t-th approximation
        row_id = J * (t - 1) + 1
        u_df = u_df_full[row_id:(row_id+J-1),]
        
        ## (2) fit the regression tree via rpart()
        u_rpart = rpart(psi_u ~ ., u_df)
        
        ## (3) process the fitted tree
        
        # (3.1) obtain the (data-defined) support for each of the parameters
        
        param_support = extractSupport(u_df, D) # hybrid_approx_v1.R
        
        # (3.2) obtain the partition
        u_partition = extractPartition(u_rpart, param_support) 
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star_min(u_rpart, u_df, u_partition, D) # partition.R
        
        # ----------------------------------------------------------------------
        n_partitions = nrow(u_partition) # number of partitions 
        
        # ----------------------------------------------------------------------
        
        K = nrow(u_partition)
        
        # new declarations here: additional storage for all 3 approximations
        const_approx   = numeric(K)       # store approx that uses 1-term taylor
        # taylor_approx  = numeric(K)     # store approx that uses 2-term taylor
        # hybrid_approx  = numeric(K)     # store approx that uses both
        
        # declare terms that will be used in the log-sum-exp trick
        eta_k = numeric(K) # log of the area of each partition A_k
        
        ck_1 = numeric(K)
        
        # ----------------------------------------------------------------------
        
        # (4) compute closed form integral over each partition
        for (k in 1:n_partitions) {
            
            # print(k)
            
            # extract "representative point" of the k-th partition
            # u = param_out[k, star_ind] %>% unlist %>% unname
            
            # compute the following for log-sum-exp trick
            # ck_1[k] = -psi(u, prior)
            
            ck_1[k] = -param_out$psi_star[k]
            
            # ck_2[k] = sum(l_k * u)
            # ck_3 = numeric(D)
            
            # ------------------------------------------------------------------
            
            for (d in 1:D) {
                
                # find column id of the first lower bound
                col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
                col_id_ub = col_id_lb + 1
                
                # limits of integration, length of the interval for param d
                upper = param_out[k, col_id_ub]
                lower = param_out[k, col_id_lb]
                
                eta_k[k] = eta_k[k] + log(upper - lower)
                
            } # end of loop computing each of 1-dim integrals
            
            const_approx[k]  = ck_1[k] + eta_k[k]
            # taylor_approx[k] = ck_1[k] + ck_2[k] + sum(ck_3)
            
        } # end of for loop over the K partitions
        
        # update approximations
        const_vec[t]  = log_sum_exp(const_approx) 
        
        
        psi_star_df = param_out %>% dplyr::select(leaf_id, psi_star)
        
        u_df = u_df %>% mutate(leaf_id = u_rpart$where)
        u_df = merge(u_df, psi_star_df, by = 'leaf_id')
        
    } # end of N_approx outer loop
    
    
    return(list(const_vec    = const_vec, 
                const_approx = const_approx,   # used to determine logML approx
                n_partitions = n_partitions,
                u_df_fit     = u_df,
                param_out    = param_out,
                u_rpart      = u_rpart))
    
} 
# end of hml_const() function --------------------------------------------------








