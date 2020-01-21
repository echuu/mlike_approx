


log_det = function(xmat) {
    return(c(determinant(xmat, logarithm = T)$modulus))
}



preprocess = function(post_samps, D, prior) {
    
    
    psi_u = apply(post_samps, 1, psi, prior = prior) %>% unname() # (J x 1)
    
    # (1.2) construct u_df -- this will require some automation for colnames
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(post_samps, psi_u) # J x (D + 1)
    names(u_df) = u_df_names
    
    
    return(u_df)
    
}



approx_lil = function(N_approx, D, u_df_full, J, prior) {
    
    def_approx = numeric(N_approx)
    
    for (t in 1:N_approx) {
        
        # if (t %% 10 == 0) {
        #     print(paste("iter", t))
        # }
        
        row_id = J * (t - 1) + 1
        
        u_df = u_df_full[row_id:(row_id+J-1),]
        
        ## (2) fit the regression tree via rpart()
        u_rpart = rpart(psi_u ~ ., u_df)
        
        ## (3) process the fitted tree
        
        # (3.1) obtain the (data-defined) support for each of the parameters
        param_support = matrix(NA, D, 2) # store the parameter supports row-wise
        
        for (d in 1:D) {
            param_d_min = min(u_df[,d])
            param_d_max = max(u_df[,d])
            
            param_support[d,] = c(param_d_min, param_d_max)
        }
        
        # paste code back here
        # (3.2) obtain the partition --- moment of truth!!
        u_partition = paramPartition(u_rpart, param_support)  # partition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition, D)
        
        n_partitions = nrow(u_partition)
        c_k = numeric(n_partitions)
        zhat = numeric(n_partitions)
        
        for (k in 1:n_partitions) {
            
            star_ind = grep("_star", names(param_out))
            u = param_out[k, star_ind] %>% unlist %>% unname
            
            c_k[k] = exp(-psi(u, prior)) # (1 x 1)
            
            l_k = lambda(u, prior)
            
            integral_d = numeric(D) # store each component of the D-dim integral 
            
            for (d in 1:D) {
                
                # updated 1/14: find column id of the first lower bound
                col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
                col_id_ub = col_id_lb + 1
                
                # d-th integral computed in closed form
                integral_d[d] = - 1 / l_k[d] * 
                    exp(- l_k[d] * (param_out[k, col_id_ub] - param_out[k, col_id_lb]))        
                
            }
            
            zhat[k] = prod(c_k[k], integral_d)
        }
        
        # return(zhat)
        
        def_approx[t] = log(sum(zhat))
        
        if (is.nan(def_approx[t])) {
            def_approx[t] = log(-sum(zhat))
        }
        
        
        
    }
    
    return(def_approx)
    
    # return(0)
    
} # end of approx_lil()

