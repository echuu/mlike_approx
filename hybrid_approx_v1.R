


hybrid_mlik = function(N_approx, D, u_df_full, J, prior) {
    
    
    hybrid_vec = numeric(N_approx)
    taylor_vec = numeric(N_approx)
    
    # compute approximation to LIL N_approx times
    for (t in 1:N_approx) {
        
        # if (t %% 10 == 0) {
        #     print(paste("iter", t))
        # }
        
        ## (1) subset out rows in u_df_full to be used in the t-th approximation
        row_id = J * (t - 1) + 1
        u_df = u_df_full[row_id:(row_id+J-1),]
        
        ## (2) fit the regression tree via rpart()
        u_rpart = rpart(psi_u ~ ., u_df)
     
        # (3.1) obtain the (data-defined) support for each of the parameters
        param_support = extractSupport(u_df, D)
        
        # (3.2) obtain the partition
        u_partition = extractPartition(u_rpart, param_support)  # partition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition, D)
        
        n_partitions = nrow(u_partition)     # numebr of partitions 
        c_k          = numeric(n_partitions) # constant term for k-th partition
        zhat         = numeric(n_partitions) # integral over k-th partition
        
        
        K = nrow(u_partition) # number of partitions 
        
        # newly declared storage
        partition_integral = numeric(K) # store the numerical integral for each A_k
        area_k             = numeric(K) # store the area of each partition
        
        taylor1_approx = numeric(K)     # store approx that uses 1-term taylor
        taylor2_approx = numeric(K)     # store approx that uses 2-term taylor
        hybrid_approx  = numeric(K)     # store approx that uses both above approx
        
        e_ck_1   = numeric(K)           # store first constant in taylor approx
        e_ck_2   = numeric(K)           # store second constant in taylor approx
        
        
        lambda_k = data.frame(matrix(NA, K, D)) # store gradient at u_k_star
        names(lambda_k) = c("lambda1", "lambda2")
        
        taylor2_integral = numeric(K)   # store integral of e^(-lambda_k'u) over A_k
        
        
        star_ind = grep("_star", names(param_out)) # columns of u_k_star for each k
        
        
        for (k in 1:K) {
            
            u1_b = param_out[k, 6]  # upper bound of u1
            u1_a = param_out[k, 5]  # lower bound of u1
            
            u2_b = param_out[k, 8]  # upper bound of u2
            u2_a = param_out[k, 7]  # lower bound of u2
            
            # (0) true value of the integral over the k-th partition
            # result = integral2(psi_fun, xmin = u1_a, xmax = u1_b, 
            #                   ymin = u2_a, ymax = u2_b, reltol = 1e-50)
            # partition_integral[k] = result$Q
            
            # ----------------------------------------------------------------------
            
            # for each parttion, A_k, keep track of the following quantities: 
            #     (1) one term taylor approximation
            #         (1.1) e_ck_1 = e^(-psi(u_kstar))
            #         (1.2) simple approx: e_ck_0 * (Area of A_k partition)
            #     (2) two term taylor approximation
            #         (2.1) lambda_k = gradient evaluated at u_k_star
            #         (2.2) e_ck_2 = e^(lambda_k'u_k_star) = 2nd (constant) term
            #         (2.3) taylor2_integral = integral of e^(-lambda_k'u) over A_k
            #         (2.4) taylor2_approx = e_ck_1 * e_ck_2 * taylor2_integral
            
            # (1) compute integral via one term taylor approximation
            u_k_star = param_out[k, star_ind] %>% unlist %>% unname
            
            # e_ck[k] = exp(-psi(u_k_star, prior))
            area_k[k] = (u1_b - u1_a) * (u2_b - u2_a)
            
            e_ck_1[k] = exp(-psi(u_k_star, prior))
            
            taylor1_approx[k] = e_ck_1[k] * area_k[k]
            
            # (2) compute integral via two term taylor approximation
            
            lambda_k[k,] = lambda(u_k_star, prior)
            lambda1  = lambda_k[k, 1]
            lambda2  = lambda_k[k, 2]
            
            
            # (2.2) e_ck_2 = e^(lambda_k'u_k_star) = 2nd (constant) term
            e_ck_2[k] = exp(sum(lambda_k[k,] * u_k_star))
            
            # (2.3) taylor2_integral = integral of e^(-lambda_k'u) over A_k
            # normally this is D-iteration loop that calculates each 1-d integral
            taylor2_integral[k]  = - 1 / lambda1 * - 1 / lambda2 * 
                (exp(-lambda1 * u1_b) - exp(-lambda1 * u1_a)) * 
                (exp(-lambda2 * u2_b) - exp(-lambda2 * u2_a))
            
            taylor2_approx[k] = e_ck_1[k] * e_ck_2[k] * taylor2_integral[k]
            
        }
        
        # take out the next chunk when running the algo for asymptotics --------
        
        # logZ_numer   = log(sum(partition_integral))
        logZ_taylor1 = log(sum(taylor1_approx))
        logZ_taylor2 = log(sum(taylor2_approx))
        
        all_integrals = cbind(numer = NA, 
                              taylor1 = taylor1_approx) %>% 
            cbind(taylor2 = taylor2_approx)
        
        diagnostics = all_integrals %>% cbind(lambda_k) %>% cbind(e_ck_2) %>% 
            cbind(taylor2_int = taylor2_integral)
        
        # partition_info = (param_out %>% mutate(perc_mem = n_obs / sum(n_obs))) %>% 
        #     dplyr::select(perc_mem, psi_hat, u1_star, u2_star) %>% 
        #     cbind(diagnostics) %>% arrange(desc(perc_mem)) %>% round(4)
        # 
        # partition_info = (param_out %>% mutate(perc_mem = n_obs / sum(n_obs))) %>% 
        #     dplyr::select(perc_mem, psi_hat) %>% 
        #     cbind(diagnostics) %>% arrange(desc(perc_mem))
        
        verbose_partition = (param_out %>% mutate(perc_mem = n_obs / sum(n_obs))) %>% 
            cbind(diagnostics) %>% arrange(desc(perc_mem))
        
        # ----------------------------------------------------------------------
        
        
        u_df = u_df %>% mutate(leaf_id = u_rpart$where, 
                               const_approx = 0,  const_resid = 0, 
                               taylor_approx = 0, taylor_resid = 0)
        
        partition_id = u_rpart$where %>% unique
        n_partitions = length(table(u_df$leaf_id))
        
        for (i in 1:K) {
            
            k = partition_id[i]
            
            u_k_star = param_out %>% filter(leaf_id == k) %>% select(star_ind) %>% 
                unname %>% unlist
            
            # compute constant approximation for psi
            u_df[u_df$leaf_id == k,]$const_approx = psi(u_k_star, prior) %>% c()
            
            # compute order 1 taylor approximation for psi
            # 
            diff_k = sweep(u_df %>% filter(leaf_id == k) %>% select(u1, u2), 2, 
                           FUN = '+', -u_k_star)             
            
            u_df[u_df$leaf_id == k,]$taylor_approx = c(psi(u_k_star, prior)) + 
                as.matrix(diff_k) %*% lambda(u_k_star, prior)
            
            u_df = u_df %>% mutate(const_resid  = psi_u - const_approx,
                                   taylor_resid = psi_u - taylor_approx)
        }
        
        error_df = data.frame(leaf_id = partition_id, 
                              const_sq_error = 0, taylor_sq_error = 0)
        
        for (i in 1:n_partitions) {
            
            k = partition_id[i]
            
            # const_sq_error
            sse_const  = sum(u_df[u_df$leaf_id == k,]$const_resid^2)
            sse_taylor = sum(u_df[u_df$leaf_id == k,]$taylor_resid^2)
            
            # for each partition, calulcate sum of residuals for const, taylor
            error_df = error_df %>% 
                mutate(const_sq_error = replace(const_sq_error, partition_id == k, 
                                                sse_const),
                       taylor_sq_error = replace(taylor_sq_error, partition_id == k,
                                                 sse_taylor))
        }
        
        # visualize the approximations side by side with associated SSE
        partition_approx = verbose_partition %>% 
            select(leaf_id, numer, taylor1, taylor2)
        
        partition_approx = merge(partition_approx, error_df, by = "leaf_id")
        
        # extract leaf id for which we use the taylor approximation
        taylor_index = error_df %>% filter(taylor_sq_error < const_sq_error) %>% 
            select(leaf_id) %>% unlist %>% unname
        
        # extract leaf id for which we use the constant approximation
        const_index = error_df %>% filter(taylor_sq_error >= const_sq_error) %>% 
            select(leaf_id) %>% unlist %>% unname
        
        
        # partition_info %>% select(numer) %>% sum %>% log
        
        const_contribution = verbose_partition %>% 
            filter(leaf_id %in% const_index) %>% 
            select(taylor1)
        
        taylor_contribution = verbose_partition %>% 
            filter(leaf_id %in% taylor_index) %>% 
            select(taylor2)
        
        # check if either contribution is empty
        if (dim(taylor_contribution)[1] == 0) {
            hybrid_approx = log(sum(const_contribution)) # 2.113237
        } else if (dim(const_contribution)[1] == 0) {
            hybrid_approx = log(sum(taylor_contribution)) # 2.113237
        } else {
            hybrid_approx = log(sum(const_contribution, taylor_contribution)) # 2.113237
        }
        
        hybrid_vec[t] = hybrid_approx
        taylor_vec[t] = logZ_taylor2
           
    } # end of outer for loop
    
    
    return(list(hybrid_vec = hybrid_vec, taylor_vec = taylor_vec))
    # return(approx_vec)
    
} # end hybrid_mlik() function




