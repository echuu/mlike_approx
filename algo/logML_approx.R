
# logML_approx.R

## functions in this file
##
##     hml_const()
##     hml()
##


#### hml_const() ---------------------------------------------------------------
#
#
hml_const = function(N_approx, D, u_df_full, J, prior) {
    
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
        param_support = extractSupport(u_df, D) #
        
        # (3.2) obtain the partition
        u_partition = extractPartition(u_rpart, param_support) 
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star_hat(u_rpart, u_df, u_partition, D) # partition.R
        
        # ----------------------------------------------------------------------
        n_partitions = nrow(u_partition) # number of partitions 
        
        # ----------------------------------------------------------------------
        
        K = nrow(u_partition)
        
        # new declarations here: additional storage for all 3 approximations
        const_approx   = numeric(K)       # store approx that uses 1-term taylor
        
        # declare terms that will be used in the log-sum-exp trick
        eta_k = numeric(K) # log of the area of each partition A_k
        
        ck_1 = numeric(K)
  
        # ----------------------------------------------------------------------
        
        # (4) compute closed form integral over each partition
        for (k in 1:n_partitions) {
            
            # print(k)
            
            # compute the following for log-sum-exp trick
            # ck_1[k] = -psi(u, prior)
            
            ck_1[k] = -param_out$psi_star[k]
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
    
    
    return(list(const_vec     = const_vec, 
                const_approx  = const_approx,   # used to determine logML approx
                n_partitions  = n_partitions,
                u_df_fit      = u_df,
                param_out     = param_out,
                u_rpart       = u_rpart,
                param_support = param_support))
    
} 
# end of hml_const() function --------------------------------------------------


# end of logML_approx.R
