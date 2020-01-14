

library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function
library('microbenchmark')

# path for lenovo
# setwd("C:/Users/ericc/mlike_approx")

# path for dell
setwd("C:/Users/chuu/mlike_approx")
source("partition/partition.R")      # load partition extraction functions
source("mvn_ig_helper.R")  # load functions specific to this model

set.seed(1)

D = 10           # dimension of paramter
p = D - 1        # dimension of beta
N = 50           # sample size

mu_beta = rep(0, p)      # mu_beta
V_beta = diag(1, p)      # V_beta
a_0 = 2 / 2              # a
b_0 = 1 / 2              # b 

beta    = c(5, -1)     # (p x 1) true coefficient vector

beta    = sample(-10:10, p, replace = T)

if (length(beta) != p) {
    print("dimension of beta must equal p")
}

sigmasq = 4            # (1 x 1) true variance

I_N = diag(1, N)       # (N x N) identity matrix
I_p = diag(1, p)       # (p x p) identity matrix

X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = t(rmvnorm(1, mean = rep(0, N), sigma = sigmasq * I_N)) # (N x 1)

y = X %*% beta + eps # (N x 1) response vector


# compute posterior parameters
V_beta_inv = solve(V_beta)
V_star_inv = t(X) %*% X + V_beta_inv

V_star  =  solve(V_star_inv)                                  # (p x p)
mu_star =  V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta)   # (p x 1)
a_n     =  a_0 + N / 2                                        # (1 x 1)
b_n     =  b_0 + 0.5 * (t(y) %*% y +                          # (1 x 1)
                            t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                            t(mu_star) %*% V_star_inv %*% mu_star) %>%  c()
# sanity check: posterior mean of beta, sigmasq (values make sense!)
mu_star
b_n / (a_n - 1)

# create prior, posterior objects
prior = list(V_beta = V_beta, 
             mu_beta = mu_beta, 
             a_0 = a_0, 
             b_0 = b_0,
             y = y, X = X,
             V_beta_inv = V_beta_inv)

# store posterior parameters
post  = list(V_star  =  V_star,
             mu_star =  mu_star,
             a_n     =  a_n,
             b_n     =  b_n,
             V_star_inv = V_star_inv)

# posterior = list(mu_star = mu_star, V_star = V_star, a_n = a_n, b_n = b_n)

# compute true log marginal likelihood
LIL_mvn_ig = lil(y, X, prior, post)

print(LIL_mvn_ig) # -106.3046 ///// -110.9457


## obtain hybrid marginal likelihood approximation

set.seed(1)
N_approx = 1 # number of approximations to compute
def_approx = approx_lil(N_approx, prior, post, D) # (100 x 1)

def_approx # -111.433

# evaluate the mean, variance of the approximations
mean(def_approx, na.rm = TRUE) # -106.1804
var(def_approx, na.rm = TRUE)  # 4.491305
print(LIL_mvn_ig)










# ------------------------------------------------------------------------------


source("mvn_ig_helper.R")  # load functions specific to this model
microbenchmark(
    def_approx = approx_lil(N_approx, prior, post, D), # (100 x 1)
    times = 5
)




# ------------------------------------------------------------------------------

## TODO: perform algorithm for a single approximation (see code starting at 
##       line 144 in nig_2d_vec.R for skeleton)

J = 3000

## (0) sample from posterior (assumed that we're able to do this) --------------
# set.seed(2)

# (0.1) sample from sigmasq | y
sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n) # (J x 1)

# (0.2) sample from mu | sigmasq, y
beta_post = data.frame(matrix(0, J, p))                         # (J x p)

for (j in 1:J) {
    # each posterior sample is stored row-wise
    beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
}

# store the posterior samples (row-wise) in a (J x D) matrix
# colnames: ( sigmasq_post, beta_post.X1, beta_post.X2, ... , beta_post.Xp )
# ** order of paramters matters, since helper functions assume a 
#    certain order when defined
u_post = data.frame(beta_post = beta_post, sigmasq_post = sigmasq_post)



# ------------------------------------------------------------------------------

## (1) create u_df (required by the functions in partition.R)

# START TEST : move this testing code elsewhere later 

# (u_post_test = head(u_post))

# test for vector input
# psi_true_mvn(u_post_test[2,], post = post)

# test for matrix input, using apply()
# apply(u_post_test, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)

## END TEST -------


# (1.1) compute psi_true() to be passed into the tree

psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)

# head(psi_u)


# (1.2) construct u_df -- this will require some automation for colnames
u_df_names = character(D + 1)
for (d in 1:D) {
    u_df_names[d] = paste("u", d, sep = '')
}
u_df_names[D + 1] = "psi_u"

# populate u_df
u_df = cbind(u_post, psi_u) # J x (D + 1)

# rename columns (needed since these are referenced explicitly in partition.R)
names(u_df) = u_df_names

# use for testing
# u_post_test = head(u_post)

# apply(u_post_test, 1, psi_true_mvn, post = post) %>% unname()

# log_mvnig(u_post_test[6,], post) 

# psi_true_mvn(u_post_test[6,], post) 

# ------------------------------------------------------------------------------

## (2) fit the regression tree via rpart()

u_rpart = rpart(psi_u ~ ., u_df)

# view the splits in the tree
# plot(u_rpart)

# START TEST : move this testing code elsewhere later 

# psi_mvn() used in main algo

# test for vector input
# psi_mvn(u_post_test[1,], prior)

# test for matrix input, using apply()
# apply(u_post_test, 1, psi_mvn, prior = prior) %>% unname() # (J x 1)

## END TEST -------

# ------------------------------------------------------------------------------

## (3) process the fitted tree

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition --- moment of truth!!
u_partition = paramPartition(u_rpart, param_support)  # partition.R

head(u_partition)

# (3.3) organize all data into single data frame (see partition.R for format)

# extracts u_star, representative point of each partition (u_star \in R^D)
# psi_hat, leaf_id, u1_star, u2_star, ... , uD_star, 
#                   u1_lb, u1_ub, ...uD_lb, uD_ub
param_out = u_star(u_rpart, u_df, u_partition, D)

head(param_out)

# ------------------------------------------------------------------------------

## DONE 1/14 : define, validate lambda() --> mvn_ig_helper.R
## DONE 1/14 : test/validate lambda (closed form vs. numerical) 

# testing: psi_mvn(), lambda_mvn_closed(), lambda_mvn()

k = 2
u = c(param_out[k,]$u1_star, param_out[k,]$u2_star, param_out[k,]$u3_star)

exp(-psi_mvn(u, prior))     # (1 x 1) -- c_k[k] calculation

lambda_mvn_closed(u, prior) # (D x 1) -- lambda(u_star) calculation

# check that numerical evaluation of gradient matches above closed form
lambda_mvn(u, prior) # (D x 1)

# check speed of both of these functions
#microbenchmark(
#    closed = lambda_mvn_closed(u, prior),
#    numerical = lambda_mvn(u, prior),
#    times = 50
#)

# closed form is about 120x faster



## (4) begin main algorithm 
n_partitions = nrow(u_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)

for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    # c_k[k] = exp(-psi(param_out[k,]$u1_star, 
    #                   param_out[k,]$u2_star,
    #                   y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    # TODO: avoid explicit definition of each representative point
    u = c(param_out[k,]$u1_star, param_out[k,]$u2_star, param_out[k,]$u3_star)
    
    c_k[k] = exp(-psi_mvn(u, prior)) # (1 x 1)
    
    l_k = lambda_mvn_closed(u, prior)
    
    integral_d = numeric(D) # store each component of the D-dim integral 
    
    # nothing to refactor in this loop (i think?) since we're just iterating
    # thru each of the integrals and computing an exponential term
    for (d in 1:D) {
        
        # verify these -- these need to be recalculated if the form of param_out
        # changes (if columns get shuffled)
        
        # col id will change for D > 2
        # DONE: generalize this better so there's less obscure calculation
        # col_id_lb = 5 + 2 * (d - 1)
        # col_id_ub = col_id_lb + 1
        
        # updated 1/14: find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out))
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
            exp(- l_k[d] * (param_out[k, col_id_ub] - param_out[k, col_id_lb]))        
        
    }
    
    zhat[k] = prod(c_k[k], integral_d)
}

log(sum(-zhat))
print(LIL_mvn_ig)


# end of single approximation simulation ---------------------------------------

## DONE : perform algorithm for a batch approximation (averaged) -- 
##        essentially just pasting the code above into a loop






