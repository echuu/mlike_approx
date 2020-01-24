
extractPartition = function(u_tree, param_support = NULL) {
    
    # pre-process the partition output
    rules = rpart.rules(u_tree, digits = 6, roundint = FALSE)
    
    ## 1/23 -- updated following 3 lines
    rules_df_full = apply(rules, 2, as.character)
    psi_hat_rules = round(as.numeric(rules_df_full[,1]), 4)
    rules_df = rules_df_full[,-c(1:2)]
    rules_str = data.frame(x = matrix(0, nrow(rules_df))) # n_partitions x 1
    
    # predicted values for each partition
    # TODO: generalize the 4 to user input later
    # psi_hat = round(as.numeric(rules[,1]), 4) 
    
    partition_id = sort(unique(u_tree$where)) # row id of leaf node information
    psi_hat = round(u_tree$frame[partition_id,]$yval, 4)
    
    # n_params = length(u_tree$variable.importance)
    
    n_params = length(u_tree$ordered)
    n_partitions = nrow(rules_str)
    
    ## 1/23 -- added the following check
    # should include this check to make sure the values match up
    if (sum(psi_hat_rules %in% psi_hat) != n_partitions) {
        print("warning -- one or more partition rules lost when rounding")
    }
    
    for(r in 1:nrow(rules)) {
        rules_str[r,] = str_c(rules_df[r,][rules_df[r,] != ""], collapse = ' ')
    }
    
    # use (-Inf, Inf) as the support for all the parameters if not specified
    # by user
    if (is.null(param_support)) {
        param_support = cbind(rep(-Inf, n_params), rep(Inf, n_params))
    }
    
    # initialize storage for the partition intervals -- (n_params x 2)
    partition = data.frame(matrix(NA, n_partitions, n_params * 2))
    
    for (p in 1:n_params) {
        
        p_col = 2 * p - 1
        
        # form column names for {(u1_lb, u1_ub),...,(up_lb, up_ub)}
        lb = paste("u", p, "_lb", sep = "")
        ub = paste("u", p, "_ub", sep = "")
        names(partition)[p_col:(p_col + 1)] = c(lb, ub)
        
        # fill in the support for the p-th parameter
        partition[, p_col] = param_support[p, 1]       # param p lower bound
        partition[, p_col + 1] = param_support[p, 2]   # param p upper bound
    } 
    
    
    # populate the partition with decision rules
    for (r in 1:n_partitions) {
        
        # (0)   split the string by & (each decision stored in matrix/vector)
        part_vec = str_split(rules_str[r,], "\\&", simplify = TRUE)
        
        # (1)   iterate over the partition/decisions
        for (p in 1:ncol(part_vec)) {
            
            
            ### fixed : node_partition() function is buggy
            processNode = node_partition(part_vec[p], param_support)
            
            # (1.1) identify the parameter corresponding to each partition
            param_id = processNode$id
            col_id = 2 * param_id - 1
            
            # (1.2) extract values defining the partition/decision
            bdry = processNode$interval
            
            # (1.3) store the values from (1.2) into correct column
            partition = updatePartition(partition, r, col_id, bdry)
            
        } # end of update step for each of the parameters
        
    } # end of loop storing the parititon boundaries
    
    ## 1/23 -- updated the left-appended column
    partition = cbind(psi_hat = psi_hat_rules, partition)
    
    # number of observations in each leaf node
    part_obs_tbl = table(u_tree$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    #### (2) obtain predicted value for each of the observations
    psi_hat_leaf = cbind(leaf_id = partition_id,
                         psi_hat = u_tree$frame[partition_id,]$yval) %>% 
        data.frame()
    
    partition_out = merge(psi_hat_leaf %>% 
                              mutate(psi_hat = round(psi_hat, 4)), 
                          partition, by = 'psi_hat')
    
    return(partition_out)
    
} # end paramPartition() function ----------------------------------------------



