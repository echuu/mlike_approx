

sup_u1 = c(-3, 10)
sup_u2 = c(-4, 5)


sup_u1 = c(-Inf, Inf)
sup_u2 = c(-Inf, Inf)

param_support = rbind(sup_u1, sup_u2)


## test the functions written in partition.R

# load the truncated bivariate normal data


library("tmvtnorm")   # truncated normal
library("data.table") # data handling
library("dplyr")      # data handling
library("ggplot2")    # visualizations
library("tree")

library('rpart')
library('rpart.plot')


set.seed(123)

mu = c(0.5, 0.5)
Sigma = matrix(c(1, 0.8, 0.8, 2), 2, 2)

# number of samples to draw
J = 10000

# sample from bivariate normal
X = rmvnorm(J, mu, Sigma)
plot(X)

# truncation points: x1 truncated to (-1, 0.5), x2 truncated to (-inf, 4)
a = c(-1, -Inf) # lower truncation
b = c(0.5, 4)   # upper truncation

# define psi() function --------------------------------------------------------
psi = function(u, mu, Sigma, a, b) {
    loglik = dtmvnorm(u, mu, Sigma, a, b, log = TRUE)
    return(loglik)
} # end psi() function


# sample J = 1000 using rejection sampling
X_r = rtmvnorm(J, mean = mu, sigma = Sigma, lower = a, upper = b, 
               algorithm = 'rejection') # 1000 x 2
plot(X_r)

# sample J = 1000 using gibbs sampling --> we use this to fit the tree
X_g = rtmvnorm(J, mean = mu, sigma = Sigma, lower = a, upper = b, 
               algorithm = 'gibbs')     # 1000 x 2
plot(X_g)


# evaluate psi(u)
psi_u = psi(X_g, mu, Sigma, a, b)                            # (J x 1)
u_df = data.frame(u1 = X_g[,1], u2 = X_g[,2], psi_u = psi_u) # (J x 3)


# fit decision tree using tree()
u_tree = tree(psi_u ~ u1 + u2, u_df)
plot(u_tree)
text(u_tree, cex = 0.5)

# fit decision tree using rpart()
u_rpart = rpart(psi_u ~ u1 + u2, u_df)
plot(u_rpart)
text(u_rpart, cex = 0.8)

# tree() and rpart() fits are very similar just by looking at the plotted 
# partition trees

# output of the partitions for the rpart tree
tbiv_rules = rpart.rules(u_rpart)

# obtain partition so that we can use it more conveniently
tbiv_supp = cbind(a, b) # u1 truncated to (-1, 0.5), u2 truncated to (-Inf, 4)
paramPartition(u_rpart, tbiv_supp)

# -----------------------------------------------------------------------------

# TODO: find the decision rules for each of the observations (params)

# TODO: for subsets of the observations (params), obtain another tree (and 
# thus a finer partition for these observations)



# in progress ------------------------------------------------------------------

# rationale: we calculate the value of the objective function for each point 
# inside a partition (this is actually already calculated) and declare the 
# point to be representative of the partition whose *objective function value 
# is the median objective function value in that partition.*

# given fitted rpart() object, obtain evaluated objective function for each
# of the observations

# see link below to see descriptions of the return of rpart object
# https://stat.ethz.ch/R-manual/R-devel/library/rpart/html/rpart.object.html

# where: an integer vector of the same length as the number of observations 
# in the root node, containing the row number of frame corresponding to the leaf 
# node that each observation falls into. 

# idea: compute the predicted value for each of values (this will result in 
# same values for points in the same partition/rectangle), subtract from true
# value (psi(u), for each u in R^2), square -> objective function value

# checks: the sum of these residuals for EACH partition should match deviation 
# provided in the rpart.object$frame output in the 'dev' column

# steps to extract the representative points of each partition

#### (1) obtain partition location of each observation

# (1.1) determine which partition each observation is grouped in
u_df = u_df %>% mutate(leaf_id = u_rpart$where)

# (1.2) obtain the rows in rpart.object$frame associated with the leaf nodes
# these rows contain the fitted value for nodes that fall w/in a given partition
leaf_id = sort(unique(u_rpart$where)) # row id of leaf node information

# number of observations in each leaf node
part_obs_tbl = table(u_rpart$where) %>% data.frame
names(part_obs_tbl) = c("leaf_id", "n_obs")

#### (2) obtain predicted value for each of the observations
psi_hat = cbind(leaf_id, psi_hat = u_rpart$frame[leaf_id,]$yval)
u_df = merge(u_df, psi_hat, "leaf_id")

#### (3) compute squared residuals for each of the observations
u_df = u_df %>% mutate(dev_sq = (psi_u - psi_hat)^2)

#### (4) for each partition: compute the median of the squared residuals 

# check: residual for each of the nodes should match the deviance in 'frame'
u_df %>% group_by(leaf_id) %>% summarise(sum(dev_sq)) # matches!!!

# rownames(u_df) = NULL

## at this point u_df looks like: ----------------------------------------------
# | leaf_id |    u1   |    u2   | ... |   up  |  psi_u  |  psi_hat  |  dev_sq
# |---------------------------------------------------------------------------
# |       4 |  0.423  | -4.584  | ... |   up  | -10.436  |  -6.522  |  15.315
# |       4 | -0.425  | -4.455  | ... |   up  | -8.1148  |  -6.522  |  2.5353
# |     ... |    ...  |    ...  | ... |   up  |    ...   |     ...  |    ... 
## ----------------------------------------------------------------------------- 

# (4.1) for each partition: sort the rows by squared residual values
n_partitions = length(leaf_id)
n_params = 2 # TODO: this should be passed in at the beginning of the function

# initialize data.frame to store the representative points of each partition
# (n_partitions x n_params) 
u_star = matrix(NA, n_partitions, n_params) %>% data.frame() 

for (k in 1:n_partitions) {
    
    # number of observations in partition k
    n_obs = part_obs_tbl$n_obs[k]
    
    # subset out observations in partition k, sort on dev_sq column
    sorted_partition = u_df %>% filter(leaf_id == leaf_id[k]) %>% 
        arrange(dev_sq)
    
    # (4.2) for each partition: save the row whose squared residual value is 
    # the median; this point is that partitions's "representative point"
    
    # extract row corresponding to the median
    med_row = floor(n_obs / 2)
    part_k_med = sorted_partition[med_row,]
    
    # extract the 'representative point' of partition k -- this will be a 
    # p-dim vector
    u_vec_k = (part_k_med %>% select(u1:psi_u))[, -(n_params + 1)] %>% 
        as.matrix() %>% c()
    
    u_star[k,] = u_vec_k
} # end of loop extracting representative points



# in progress ------------------------------------------------------------------



# TODO: adjust the hyperparameters in rpart() to encourage a deeper tree (finer
# partitions) to be grown, one that is more flexible - in this case we aren't
# too worried about partitions, since we just need many partitions so that we 
# have enough points used in the approximation of the marginal likelihood

# see link below to adjust the rpart() fit to allow for deeper tree growth
# https://www.rdocumentation.org/packages/rpart/versions/4.1-15/topics/rpart.control

# TODO: consider KS distance as a means of thresholding for deciding whether or
# not we need to a finer partition

# GOAL -------------------------------------------------------------------------
# TODO: compute the marginal likelihood using MCMC techniques

# TODO: compute the (approximate) marginal likelihood using the points
# that the partitions give us in the previous steps





