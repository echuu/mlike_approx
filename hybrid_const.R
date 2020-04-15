

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
        
        param_support = extractSupport(u_df, D) # hybrid_approx_v1.R
        
        # (3.2) obtain the partition
        u_partition = extractPartition(u_rpart, param_support)  # extractPartition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R
        
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
        
        # star_ind will be a vector of indices -- subsetting these out of 
        # param_out will give u_k = (u_k1, u_k2, ... , u_kD)
        star_ind = grep("_star", names(param_out)) 
        
        # ----------------------------------------------------------------------
        
        # (4) compute closed form integral over each partition
        for (k in 1:n_partitions) {
            
            # print(k)
            
            # extract "representative point" of the k-th partition
            u = param_out[k, star_ind] %>% unlist %>% unname
            
            # compute the following for log-sum-exp trick
            ck_1[k] = -psi(u, prior)
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
                
                # ck_3[d] = log_int_rect(l_k[d], lower, upper)
                
            } # end of loop computing each of 1-dim integrals
            
            const_approx[k]  = ck_1[k] + eta_k[k]
            # taylor_approx[k] = ck_1[k] + ck_2[k] + sum(ck_3)
            
        } # end of for loop over the K partitions
        
        # update approximations
        const_vec[t]  = log_sum_exp(const_approx) 
        
        u_df = u_df %>% mutate(leaf_id = u_rpart$where,
                               const_approx = 0,  const_resid = 0)

        partition_id = u_rpart$where %>% unique
       
        # for each partition, compute the sum of squared residuals,
        # (psi(u) - psi_tilde(u))^2
        for (j in 1:K) {

            k = partition_id[j]

            u_k_star = param_out %>% filter(leaf_id == k) %>%
                dplyr::select(star_ind) %>% unname %>% unlist

            #### compute constant approximation for psi
            # note: we actually already do this when computing e_ck_1, so
            # eventually, we can just take log(e_ck_1) to recover this term
            u_df[u_df$leaf_id == k,]$const_approx = psi(u_k_star, prior) %>% c()

            #### compute order 1 taylor approximation for psi

            # assumes u1,...,uD are the first D columns of u_df -- make sure
            # this structure is maintained, maybe provide a check ?
            # diff_k = sweep(u_df %>% filter(leaf_id == k) %>%
            #                    dplyr::select(c(1:D)), 2,
            #                FUN = '+', -u_k_star)

            # methodology notes:
            # compute psi(u_star) for each of the partitions; we do this in 2
            # ways -> (1) constant approximation, (2) order 1 taylor
            # based on the residual for each approximation, we decide
            # which approximation to use to compute the integral over each
            # partition
            # u_df[u_df$leaf_id == k,]$taylor_approx = c(psi(u_k_star, prior)) +
            #     as.matrix(diff_k) %*% lambda(u_k_star, prior)

            # compute difference between psi_u and corresponding approximation
            u_df = u_df %>% mutate(const_resid  = psi_u - const_approx)

        } # end of for loop computing residuals
    
    } # end of N_approx outer loop

        
    return(list(const_vec    = const_vec, 
                const_approx = const_approx,
                n_partitions = n_partitions,
                u_df_fit     = u_df,
                param_out    = param_out))

} # end of hml_const() function





















