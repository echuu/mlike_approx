

## ** will be moving algorithm functions into this file as they are 
## they are updated to improve fit



## TODO: write function that generates an output like the following
## see mvn_troubleshoot.R for inspiration of what the output should look like

#    psi_hat perc_mem partition_integral taylor1_integral taylor2_numer
# 1   0.2875   0.4438        2.869490199       3.02938082    5.29740423
# 2   1.1606   0.1032        0.683710569       0.79907025    1.39124488
# 3   0.8198   0.0878        0.497127890       0.52349531    0.82736043
# 4   1.0690   0.0796        0.490854427       0.53376998    1.23866939
# 5   1.1197   0.0732        0.451552152       0.49167407    0.65655063
# 6   1.8453   0.0564        0.338489901       0.41924842    0.75335918
# 7   2.3347   0.0388        0.255096119       0.45286334   14.26437054
# 8   2.4110   0.0370        0.234916417       0.30812222    0.36149494
# 9   2.2413   0.0362        0.209963191       0.46940746    0.67163583
# 10  3.6971   0.0188        0.087049219       0.24374799    1.29709727
# 11  3.4981   0.0120        0.078856266       0.17221502    0.40278265
# 12  4.8799   0.0054        0.047904170       0.13637846    2.00592708
# 13  4.6801   0.0048        0.022964372       0.04403228    0.03690540
# 14  4.3749   0.0030        0.009180518       0.02178593    0.01005086





# using median objective function does worse in the MVN example 
# TODO: check performance in singular example

u_star_med = function(rpart_obj, u_df_in, partition, n_params) {
    
    # notes: 
    # previously, we considered the u_k_star to be the point whose 
    # squared residual value is the median in the k-th partition
    
    # here, we take u_k_star to be the point whose evaluation of psi() is 
    # the median value in the k-th partition (they should be the same points?)
    
    # (1.1) determine which partition each observation is grouped in
    u_df = u_df_in %>% mutate(leaf_id = rpart_obj$where)
    
    # (1.2) obtain the rows in rpart.object$frame associated with the leaf nodes
    # these rows contain the fitted value for nodes that fall w/in a 
    # given partition
    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node
    
    # number of observations in each leaf node
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    #### (2) obtain predicted value for each of the observations
    psi_hat_leaf = cbind(leaf_id = partition_id,
                         psi_hat = rpart_obj$frame[partition_id,]$yval) %>% 
        data.frame()
    
    # append the fitted value for each row on the right as an additional column
    u_df = merge(u_df, psi_hat_leaf, "leaf_id")
    
    #### (3) compute squared residuals for each of the observations
    # u_df = u_df %>% mutate(dev_sq = (psi_u - psi_hat)^2)
    
    #### (4) for each partition: compute the median of the squared residuals 
    
    # check: residual for each of the nodes should match the deviance in 'frame'
    # u_df %>% group_by(leaf_id) %>% summarise(sum(dev_sq)) # matches!!!
    
    ## at this point u_df looks like: ------------------------------------------
    #
    # | leaf_id |   u1   |   u2   | ... |   up  |  psi_u  |  psi_hat |  dev_sq
    # |------------------------------------------------------------------------
    # |       4 |  0.423 | -4.584 | ... |   up  | -10.436 |  -6.522  |  15.315
    # |       4 | -0.425 | -4.455 | ... |   up  | -8.1148 |  -6.522  |  2.5353
    # |     ... |    ... |    ... | ... |   up  |    ...  |     ...  |    ... 
    #
    ## -------------------------------------------------------------------------
    
    # (4.1) for each partition: sort the rows by squared residual values
    n_partitions = length(partition_id)
    # n_params = 2 # DONE: this should be passed into the function
    
    # initialize data.frame to store the representative points of each partition
    # (n_partitions x n_params) 
    # DONE: column names should match the convention of 'u1, u2, ... '
    # extra column for the lead id
    u_star_df = matrix(NA, n_partitions, n_params + 1) %>% data.frame() 
    
    for (k in 1:n_partitions) {
        
        # number of observations in partition k
        n_obs = part_obs_tbl$n_obs[k]
        
        # 1/24 -- CHANGED WHAT WE SORT ON
        # subset out observations in partition k, sort on psi_u column
        sorted_partition = u_df %>% filter(leaf_id == partition_id[k]) %>% 
            arrange(psi_u) 
        
        # (4.2) for each partition: save the row whose PSI value is 
        # the median; this point is that partitions's "representative point"
        
        # extract row corresponding to the median
        med_row = floor(n_obs / 2)
        part_k_med = sorted_partition[med_row,]
        
        # extract the 'representative point' of partition k -- this will be a 
        # p-dim vector
        u_vec_k = (part_k_med %>% 
                       dplyr::select(u1:psi_u))[, -(n_params + 1)] %>% 
            as.matrix() %>% c()
        
        u_star_df[k,] = c(u_vec_k, part_k_med$leaf_id)
    } # end of loop extracting representative points
    
    # 1/14 -- generalizing this to D many parameters -- DONE
    u_star_names = character(n_params + 1)
    for (d in 1:n_params) {
        u_star_names[d] = paste(names(u_df_in)[d], '_star', sep = '')
    }
    u_star_names[n_params + 1] = "leaf_id"
    
    names(u_star_df) = u_star_names
    
    # names(u_star_df) = c("u1_star", "u2_star", "leaf_id")
    
    ## merge with the boundary of each of the partitions
    u_df_full = merge(u_star_df, partition, by = 'leaf_id')
    
    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    u_df_full = merge(u_df_full, part_obs_tbl, by = 'leaf_id')
    
    
    
    
    
}



extractPartition = function(u_tree, param_support = NULL) {
    
    # pre-process the partition output
    rules = rpart.rules(u_tree, digits = 6, roundint = FALSE)
    
    ## 1/23 -- updated following 3 lines
    rules_df_full = apply(rules, 2, as.character)
    # psi_hat_rules = round(as.numeric(rules_df_full[,1]), 4)
    rules_df = rules_df_full[,-c(1:2)]
    rules_str = data.frame(x = matrix(0, nrow(rules_df))) # n_partitions x 1
    
    partition_id = sort(unique(u_tree$where)) # row id of leaf node information
    
    # this is now in the order of the rules_df_full, rules_df dataframes
    psi_hat_id = cbind(leaf_id = partition_id, 
                       psi_hat = u_tree$frame[partition_id, 5]) %>% 
        data.frame() %>% arrange(psi_hat)
    
    # psi_hat = round(u_tree$frame[partition_id,]$yval, 6)
    # psi_hat_sort = sort(psi_hat) # sort these to match the order of the rules
    
    n_params = length(u_tree$ordered)
    n_partitions = nrow(rules_str)
    
    
    for(r in 1:nrow(rules)) {
        rules_str[r,] = str_c(rules_df[r,][rules_df[r,] != ""], collapse = ' ')
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
    
    ## 1/23 -- updated the lef       t-appended column
    partition_out = cbind(psi_hat_id, partition)
    
    # number of observations in each leaf node
    # part_obs_tbl = table(u_tree$where) %>% data.frame
    # names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    #### (2) obtain predicted value for each of the observations
    # psi_hat_leaf = cbind(leaf_id = partition_id,
    #                     psi_hat = u_tree$frame[partition_id,]$yval) %>% 
    #    data.frame()
    
    #partition_out = merge(psi_hat_leaf %>% 
    #                          mutate(psi_hat = round(psi_hat, 4)), 
    #                      partition, by = 'psi_hat')
    
    return(partition_out)
    
} # end paramPartition() function ----------------------------------------------



