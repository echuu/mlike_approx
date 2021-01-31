


l1_norm = function(u, u_0) {
    sum(abs(u - u_0))
}

hyb_numer = function(u_df, psi, params,
                     D = ncol(u_df) - 1) {
    
    ## (2) fit the regression tree via rpart()
    u_rpart = rpart(psi_u ~ ., u_df)
    
    ## (3) process the fitted tree
    
    # (3.1) obtain the (data-defined) support for each of the parameters
    param_support = extractSupport(u_df, D) #
    
    # (3.2) obtain the partition
    u_partition = extractPartition(u_rpart, param_support) 
    
    
    # param_out = u_star_cand(u_rpart, u_df, u_partition, D) # partition.R
    param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R
    # opt_part = param_out$optimal_part
    
    # ----------------------------------------------------------------------
    n_partitions = nrow(u_partition) # number of partitions 
    
    # ----------------------------------------------------------------------
    
    psi_partition = param_out %>% 
        dplyr::select(-c('leaf_id', 'psi_choice', 'logQ_cstar', 'n_obs'))
    
    bounds = psi_partition %>% dplyr::select(-c("psi_star"))
    
    # --------------------------------------------------------------------------
    
    #### extension starts here -------------------------------------------------
    
    
    ### (1) find global mean
    u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean
    
    
    ### (2) find point in each partition closest to global mean (for now)
    
    # u_k for each partition
    u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)
    
    l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
    u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)
    
    # take min result, group_by() leaf_id
    psi_df = u_df_part %>% 
        group_by(leaf_id) %>% filter(l1_cost == min(l1_cost)) %>% 
        data.frame
    
    
    K = nrow(bounds)
    log_terms = numeric(K) # store terms so that we can use log-sum-exp()
    G_k = numeric(K)       # store terms coming from gaussian integral
    
    
    ##
    ### (4) compute hessian, H_k for each of the u_k
    # H_k = Q_beta # same hessian for all partitions A_k
    
    # H_k_inv = Q_beta_inv
    
    ### (5) compute mean component for the gaussian, m_k for each partition
    # m_k = mu_beta # same for all partitions A_k
    
    # lambda_k = apply(psi_df[,1:D], 1, lambda, prior = prior) # column-wise, (D x K)
    ##
    
    for (k in 1:K) {
        
        u_k = unname(unlist(psi_df[k,1:D]))
        # diff_k = u_k - m_k
        
        H_k = pracma::hessian(psi, u_k, prior = params)
        H_k_inv = chol2inv(chol(H_k))
        
        lambda_k = pracma::grad(psi, u_k, prior = params)
        b_k = H_k %*% u_k - lambda_k
        m_k = H_k_inv %*% b_k
        
        lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
        ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
        G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
        
        # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])
        
        log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
            psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
            0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
        
    }
    
    log_sum_exp(log_terms)
    
    
}


