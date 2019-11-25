
# partition.R ------------------------------------------------------------------
# 


library('mvtnorm')      # multivariate normal density
library('MASS')         # mvnorm()
library('ggplot2')      # don't think we use this in these functions
library('rpart')        # rpart() for fitting the tree
library('rpart.plot')   # do we need this?
library('tidyr')        # data mangement
library('readr')        # data management
library('stringr')      # regex functions

# partition object initialization

# when creating the partition dataframe, we need:
#    (1) number of parameters
#    (2) the support of each of the parameters
#    (3) the number of partitions (extracted rpart.rules() function)

# param_support : (p x 2) with the range that each parameter stored row-wise
paramPartition = function(u_tree, param_support = NULL) {
    
    # pre-process the partition output
    rules = rpart.rules(u_tree) 
    rules_df = apply(rules, 2, as.character)[,-c(1:2)]
    rules_str = data.frame(x = matrix(0, nrow(rules_df))) # n_partitions x 1
    
    n_params = length(u_tree$variable.importance)
    n_partitions = nrow(rules_str)
    
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
    
    return(partition)
    
} # end paramPartition() function



# ------------------------------------------------------------------------------



# updatePartition() : 
# store the values from (1.2) into correct column of the final df
# partition: data frame that stores the current decisions/boundaries
# row_id   : the row corresponding to the decision to be updated
# col_id   : the TWO columns that need to be updated (parameter id)
# bdry     : the 2-d vector that has the numerical decision boundaries
updatePartition = function(partition, row_id, col_id, bdry) {
    
    partition[row_id, col_id:(col_id + 1)] = bdry 
    
    return(partition)
}


# ------------------------------------------------------------------------------

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





