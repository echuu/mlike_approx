
## outdated_star.R
## - contains outdated implementations of u_star() function
## - has come useful computatations that may be helpful later






#### u_star_min_dev() ----------------------------------------------------------
# will not work with current version of hml_const() 
# since current hml_const() does not compute using u_k_star anymore
u_star_min_dev = function(rpart_obj, u_df_in, partition, n_params) {
    
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
    u_df = u_df %>% mutate(dev_sq = (psi_u - psi_hat)^2)
    
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
        
        # subset out observations in partition k, sort on dev_sq column
        sorted_partition = u_df %>% dplyr::filter(leaf_id == partition_id[k]) %>%
            arrange(dev_sq)
        
        # sorted_partition = u_df %>% dplyr::filter(leaf_id == partition_id[k]) %>%
        #     arrange(psi_u)
        
        # (4.2) for each partition: save the row whose squared residual value is 
        # the median; this point is that partitions's "representative point"
        
        # extract row corresponding to the median
        part_k_med = sorted_partition[1,]
        
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
    
    
    return(u_df_full)
}
# end u_star() function --------------------------------------------------------




#### u_star_max() --------------------------------------------------------------
# TODO: integrate this function  with logQ() function so that the choice 
# of psi_star is the minimizer of log(Q(c_star))
u_star_max = function(rpart_obj, u_df_in, partition, n_params) {
    
    # candidates for c_star:
    # (1) constant fitted by tree
    # (2) mean of c_{kj}
    # (3) median of c_{kj}
    # (4) max of c_{kj}
    
    # (1.1) determine which partition each observation is grouped in
    u_df = u_df_in %>% mutate(leaf_id = rpart_obj$where)
    
    # (1.2) obtain the rows in rpart.object$frame associated with the leaf nodes
    # these rows contain the fitted value for nodes that fall w/in a 
    # given partition
    partition_id = sort(unique(rpart_obj$where)) # row id of leaf node
    
    # number of observations in each leaf node
    part_obs_tbl = table(rpart_obj$where) %>% data.frame
    names(part_obs_tbl) = c("leaf_id", "n_obs")
    
    ## at this point u_df looks like: ------------------------------------------
    #
    # | leaf_id |   u1   |   u2   | ... |   up  |  psi_u  |  psi_hat |  
    # |------------------------------------------------------------------------
    # |       4 |  0.423 | -4.584 | ... |   up  | -10.436 |  -6.522  |  
    # |       4 | -0.425 | -4.455 | ... |   up  | -8.1148 |  -6.522  |  2.5353
    # |     ... |    ... |    ... | ... |   up  |    ...  |     ...  |    ... 
    #
    ## -------------------------------------------------------------------------
    
    # (4.1) for each partition: sort the rows by squared residual values
    n_partitions = length(partition_id)
    
    psi_df = data.frame(leaf_id = partition_id, psi_star = NA)
    psi_star = numeric(n_partitions)    
    
    for (k in 1:n_partitions) {
        
        # number of observations in partition k
        n_obs = part_obs_tbl$n_obs[k]
        
        # find max value for k-th partition
        psi_star[k] = u_df %>% dplyr::filter(leaf_id == partition_id[k]) %>% 
            dplyr::select(psi_u) %>% max
        
    } # end of loop extracting representative points
    
    psi_df$psi_star = psi_star
    
    # 1/14 -- generalizing this to D many parameters -- DONE
    # u_star_names = character(n_params + 1)
    # for (d in 1:n_params) {
    #     u_star_names[d] = paste(names(u_df_in)[d], '_star', sep = '')
    # }
    # u_star_names[n_params + 1] = "leaf_id"
    # 
    # names(u_star_df) = u_star_names
    
    # names(u_star_df) = c("u1_star", "u2_star", "leaf_id")
    
    ## merge with the boundary of each of the partitions
    u_df_full = merge(psi_df, partition, by = 'leaf_id')
    
    # append the number of observations for each leaf node to the right
    # this is later used to determine the type of approximation to use
    u_df_full = merge(u_df_full, part_obs_tbl, by = 'leaf_id')
    
    
    return(u_df_full)
    
}  
# end of u_star_max() function ----------------------------------------------------- 





# u_star_max = function(rpart_obj, u_df_in, partition, n_params) {
#     
#     # (1.1) determine which partition each observation is grouped in
#     u_df = u_df_in %>% mutate(leaf_id = rpart_obj$where)
#     
#     # (1.2) obtain the rows in rpart.object$frame associated with the leaf nodes
#     # these rows contain the fitted value for nodes that fall w/in a 
#     # given partition
#     partition_id = sort(unique(rpart_obj$where)) # row id of leaf node
#     
#     # number of observations in each leaf node
#     part_obs_tbl = table(rpart_obj$where) %>% data.frame
#     names(part_obs_tbl) = c("leaf_id", "n_obs")
#     
#     #### (2) obtain predicted value for each of the observations
#     psi_hat_leaf = cbind(leaf_id = partition_id,
#                          psi_hat = rpart_obj$frame[partition_id,]$yval) %>% 
#         data.frame()
#     
#     # append the fitted value for each row on the right as an additional column
#     u_df = merge(u_df, psi_hat_leaf, "leaf_id")
#     
#     #### (3) compute squared residuals for each of the observations
#     u_df = u_df %>% mutate(dev_sq = (psi_u - psi_hat)^2)
#     
#     #### (4) for each partition: compute the median of the squared residuals 
#     
#     # check: residual for each of the nodes should match the deviance in 'frame'
#     # u_df %>% group_by(leaf_id) %>% summarise(sum(dev_sq)) # matches!!!
#     
#     ## at this point u_df looks like: ------------------------------------------
#     #
#     # | leaf_id |   u1   |   u2   | ... |   up  |  psi_u  |  psi_hat |  dev_sq
#     # |------------------------------------------------------------------------
#     # |       4 |  0.423 | -4.584 | ... |   up  | -10.436 |  -6.522  |  15.315
#     # |       4 | -0.425 | -4.455 | ... |   up  | -8.1148 |  -6.522  |  2.5353
#     # |     ... |    ... |    ... | ... |   up  |    ...  |     ...  |    ... 
#     #
#     ## -------------------------------------------------------------------------
#     
#     # (4.1) for each partition: sort the rows by squared residual values
#     n_partitions = length(partition_id)
#     # n_params = 2 # DONE: this should be passed into the function
#     
#     # initialize data.frame to store the representative points of each partition
#     # (n_partitions x n_params) 
#     # DONE: column names should match the convention of 'u1, u2, ... '
#     # extra column for the lead id
#     u_star_df = matrix(NA, n_partitions, n_params + 1) %>% data.frame() 
#     
#     for (k in 1:n_partitions) {
#         
#         # number of observations in partition k
#         n_obs = part_obs_tbl$n_obs[k]
#         
#         # subset out observations in partition k, sort on dev_sq column
#         # sorted_partition = u_df %>% dplyr::filter(leaf_id == partition_id[k]) %>%
#         #     arrange(dev_sq)
#         
#         sorted_partition = u_df %>% dplyr::filter(leaf_id == partition_id[k]) %>%
#             arrange(psi_u)
#         
#         # (4.2) for each partition: save the row whose squared residual value is 
#         # the median; this point is that partitions's "representative point"
#         
#         # extract row corresponding to the median
#         u_k_row = floor(n_obs / 2)
#         u_k_row = n_obs
#         part_k_med = sorted_partition[u_k_row,]
#         # part_k_med = sorted_partition[1,]
#         
#         # extract the 'representative point' of partition k -- this will be a 
#         # p-dim vector
#         u_vec_k = (part_k_med %>% 
#                        dplyr::select(u1:psi_u))[, -(n_params + 1)] %>% 
#             as.matrix() %>% c()
#         
#         u_star_df[k,] = c(u_vec_k, part_k_med$leaf_id)
#     } # end of loop extracting representative points
#     
#     # 1/14 -- generalizing this to D many parameters -- DONE
#     u_star_names = character(n_params + 1)
#     for (d in 1:n_params) {
#         u_star_names[d] = paste(names(u_df_in)[d], '_star', sep = '')
#     }
#     u_star_names[n_params + 1] = "leaf_id"
#     
#     names(u_star_df) = u_star_names
#     
#     # names(u_star_df) = c("u1_star", "u2_star", "leaf_id")
#     
#     ## merge with the boundary of each of the partitions
#     u_df_full = merge(u_star_df, partition, by = 'leaf_id')
#     
#     # append the number of observations for each leaf node to the right
#     # this is later used to determine the type of approximation to use
#     u_df_full = merge(u_df_full, part_obs_tbl, by = 'leaf_id')
#     
#     
#     return(u_df_full)
# }
# # end of u_star_max() function--------------------------------------------------



