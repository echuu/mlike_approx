library(rpart)
library(rpart.plot)

d = 5
n = 100
X = mvtnorm::rmvnorm(100, rep(0, d), diag(1, d)) # covariates
beta = 1:d
eps = rnorm(n, 0, 3)
y = X %*% beta + eps # response

tree_model = rpart(y ~ X)

plot(tree_model)


rpart.rules(tree_model) ## partition in shitty format


param_support = extractSupport(X, d) ## get the support of the data
extractPartition(tree_model, param_support) ## partition in good format



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
    
  } # end of loop storing the partition boundaries
  
  ## 1/23 -- updated the left-appended column
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
  
} 
# end extractPartition() function ----------------------------------------------




#### updatePartition() --------------------------------------------------------- 
# store the values from (1.2) into correct column of the final df
# partition: data frame that stores the current decisions/boundaries
# row_id   : the row corresponding to the decision to be updated
# col_id   : the TWO columns that need to be updated (parameter id)
# bdry     : the 2-d vector that has the numerical decision boundaries
updatePartition = function(partition, row_id, col_id, bdry) {
  
  partition[row_id, col_id:(col_id + 1)] = bdry 
  
  return(partition)
}

# end updatePartition() function -----------------------------------------------




#### node_partition() ----------------------------------------------------------
# str_in needs only have ONE of the interval identifiers, i.e., str_split()
# already needs to have been called, splitting on "\\&"
node_partition = function(str_in, param_support) {
  
  LOWER = ">="
  UPPER = "<"
  
  # regex to subset out the numbers in the decision output
  INTERVAL = "(?>-)*[[:digit:]]+\\.*[[:digit:]]*" 
  
  # subset out the numbers that appear in the decision string output    
  decision = as.numeric(unlist(regmatches(str_in, gregexpr(INTERVAL, str_in, 
                                                           perl = TRUE))))
  
  
  if (grepl(LOWER, str_in)) {          # decision: [x, ]
    
    # decision = as.numeric(str_split(str_in, LOWER, simplify = TRUE))
    # decision = as.numeric(unlist(regmatches(str_in, 
    #                                        gregexpr(INTERVAL, str_in, 
    #                                                 perl = TRUE))))
    
    param_id = decision[1]
    
    # TODO: upper bound should be the upper bound of the support
    interval = c(decision[2], param_support[param_id, 2])   # lb interval
    
  } else if (grepl(UPPER, str_in)) {   # decision: [ ,y]
    
    # decision = as.numeric(str_split(str_in, UPPER, simplify = TRUE))
    # decision = as.numeric(unlist(regmatches(str_in, 
    #                                        gregexpr(INTERVAL, str_in, 
    #                                                 perl = TRUE))))
    
    param_id = decision[1]
    
    # TODO: lower bound should be the lower bound of the support
    interval = c(param_support[param_id, 1], decision[2])  # ub interval
    
  } else {                             # decision: [x,y]
    
    # decision = as.numeric(unlist(regmatches(str_in, 
    #                                        gregexpr(INTERVAL, str_in, 
    #                                                 perl = TRUE))))
    
    param_id = decision[1]
    u_lb = decision[2]
    u_ub = decision[3]
    
    interval = c(u_lb, u_ub)
  }
  
  return(list(id = param_id, interval = interval))
}

# end node_partition() function ------------------------------------------------






