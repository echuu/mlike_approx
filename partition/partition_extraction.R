

# partition_extraction.R

# goal: extract the partition from the fitted regression tree output
# misc: using the partitions, compute the midpoint of all the partitions


# get decision tree rule/path pattern for every row of predicted dataset
# https://tinyurl.com/t5mgkce


# what we want should be a function of this solution here..
# list of splitting conditions for each terminal node
# https://tinyurl.com/tpuvaeb


# ------------------------------------------------------------------------------

# fit ctree() on gaussian mixture to see if we're getting similar partitioning 
# scheme

library('mvtnorm')   # multivariate normal density
library('MASS')      # mvnorm()
library('ggplot2')
library('tree')

set.seed(123)

d = 2
N = 100
mu1 = c(1, 2)
mu2 = c(4, 6)
Sigma1 = d / N * diag(c(1, 1))
Sigma2 = d / N * diag(c(1, 1))

# mixture weights
pi1 = 0.2
pi2 = 0.8

J = 5000

n1 = sum(runif(J) <= pi1) # draws from 1st gaussian
n2 = J - n1               # draw sfrom 2nd gaussian


# define psi() function --------------------------------------------------------
psi = function(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2) {
    
    # log-density evaluated over the rows of u -- (J x 1)
    loglik = log(pi1 * dmvnorm(u, mu1, Sigma1) + 
                     pi2 * dmvnorm(u, mu2, Sigma2))
    
    return(loglik)
    
} # end psi() function


u = rbind(mvrnorm(n1, mu1, Sigma1), mvrnorm(n2, mu2, Sigma2)) # (J x 2)

# check distribution
u_df = data.frame(u1 = u[,1], u2 = u[,2])
# ggplot(u_df, aes(u1, u2)) + geom_point() + theme_bw()

# evaluate psi(u)
psi_u = psi(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2)           # (J x 1)
u_df = data.frame(u1 = u_df$u1, u2 = u_df$u2, psi_u = psi_u) # (J x 3)

# fit decision tree using tree() function from 'tree' package
u_tree = tree(psi_u ~ u1 + u2, u_df)
plot(u_tree)
text(u_tree, cex = 0.5)

library('rpart')
library('rpart.plot')
u_tree = rpart(psi_u ~ u1 + u2, u_df)
plot(u_tree)
text(u_tree, cex = 0.5)

head(rpart.rules(u_tree))


library('stringr')

rules = rpart.rules(u_tree)
rules_df = apply(rules, 2, as.character)[,-c(1:2)]
rules_str = data.frame(x = matrix(0, nrow(rules_df)))

for(r in 1:nrow(rules)) {
    rules_str[r,] = str_c(rules_df[r,][rules_df[r,] != ""], collapse = ' ')
}


LOWER = ">="
UPPER = "<"

test = strsplit(rules_str[1,], "\\&")[[1]]
varname = grep("u1", decision)


n_partitions = nrow(rules_str)
n_params = ncol(u_df) - 1
partition = data.frame(matrix(0, n_partitions, n_params * 2))

names(partition) = c("u1_lb", "u1_ub", "u2_lb", "u2_ub")


# create the dataframe that organizes the partition in the following form:

# | u1_lb | u1_ub | u2_lb | u2_ub | --- | up_lb | up_ub |
# |  4.36 |   -   |   -   |   -   | --- |   -   |   -   | 
# |  3.79 |  4.23 |  6.40 |   -   | --- |   -   |   -   |
# |  3.79 |  4.23 |   -   |  5.70 | --- |   -   |   -   |
# |  .... |  .... |   -   |  .... | --- |  .... |  .... |
# |  3.86 |  4.15 |  5.90 |  6.10 | --- |  .... |  .... |

for (r in 1:n_partitions) {
    
    # (0)   split the string by & (each decision stored in matrix/vector)
    part_vec = str_split(test, "\\&")
    
    # (1)   iterate over the partition/decisions
    for (p in 1:ncol(part_vec)) {
        # (1.1) identify the parameter corresponding to each partition/decision
        param_id = ""
        
        # (1.2) extract values defining the partition/decision (node_partition)
        bdry = node_partition(part_vec)
        
        # (1.3) store the values from (1.2) into correct column of the final df
        partition = updatePartition(partition, param_id, bdry)
        
    } # end of update step for each of the parameters
    
    
} # end of loop storing


# updatePartition() : 
# store the values from (1.2) into correct column of the final df
# partition: data frame that stores the current decisions/boundaries
# row_id   : the row corresponding to the decision to be updated
# col_id   : the TWO columns that need to be updated (parameter id)
# bdry     : the 2-d vector that has the numerical decision boundaries
updatePartition = function(partition, row_id, col_id, bdry) {
    
    partition_next = partition[row_id, col_id:(col_id + 1)] = bdry 
    
    return(partition_next)
}



library('tidyr')
library('readr')

# str_in needs only have ONE of the interval identifiers, i.e., str_split()
# already needs to have been called, splitting on "\\&"
node_partition = function(str_in) {
    
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
        interval = c(decision[2], Inf)   # lower-bounded interval
        
    } else if (grepl(UPPER, str_in)) {   # decision: [ ,y]
        
        # decision = as.numeric(str_split(str_in, UPPER, simplify = TRUE))
        # decision = as.numeric(unlist(regmatches(str_in, 
        #                                        gregexpr(INTERVAL, str_in, 
        #                                                 perl = TRUE))))
        
        param_id = decision[1]
        interval = c(-Inf, decision[2])  # upper-bounded interval
        
    } else {                             # decision: [x,y]
        
        # decision = as.numeric(unlist(regmatches(str_in, 
        #                                        gregexpr(INTERVAL, str_in, 
        #                                                 perl = TRUE))))
        
        param_id = decision[1]
        u_lb = decision[2]
        u_ub = decision[3]
        
        interval = c(u_lb, u_ub)
    }
    
    return(list(param_id, interval))
}











