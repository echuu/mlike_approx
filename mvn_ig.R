

library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function


# setwd("C:/Users/ericc/mlike_approx/partition")

setwd("C:/Users/chuu/mlike_approx")
source("partition/partition.R")      # load partition extraction functions
source("mvn_ig_helper.R")  # load functions specific to this model

set.seed(1)

D = 3           # dimension of paramter
p = D - 1       # dimension of beta
N = 50          # sample size

mu_beta = rep(0, p)      # mu_beta
V_beta = diag(1, p)      # V_beta
a_0 = 2 / 2              # a
b_0 = 1 / 2              # b 

beta    = c(5, -1)     # (p x 1) true coefficient vector

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
prior     = list(mu_beta = mu_beta, V_beta = V_beta, a_0 = a_0, b_0 = b_0,
                 y = y, X = X)
posterior = list(mu_star = mu_star, V_star = V_star, a_n = a_n, b_n = b_n)

# compute true log marginal likelihood
LIL_mvn_ig = lil(y, X, prior, posterior)

print(LIL_mvn_ig) # -106.3046

# ------------------------------------------------------------------------------

## TODO: perform algorithm for a single approximation (see code starting at 
##       line 144 in nig_2d_vec.R for skeleton)

J = 3000

## (0) sample from posterior (assumed that we're able to do this) --------------
set.seed(1)

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


# store posterior parameters
post = list(V_star  =  V_star,
            mu_star =  mu_star,
            a_n     =  a_n,
            b_n     =  b_n)


# ------------------------------------------------------------------------------

## (1) create u_df (required by the functions in partition.R)

# START TEST : move this testing code elsewhere later 

(u_post_test = head(u_post))

# test for vector input
psi_true_mvn(u_post_test[2,], post = post)

# test for matrix input, using apply()
apply(u_post_test, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)

## END TEST -------


# (1.1) compute psi_true() to be passed into the tree

psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)




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


# START TEST : move this testing code elsewhere later 

# psi_mvn() used in main algo

# test for vector input
psi_mvn(u_post_test[1,], prior)

# test for matrix input, using apply()
apply(u_post_test, 1, psi_mvn, prior = prior) %>% unname() # (J x 1)

## END TEST -------

# ------------------------------------------------------------------------------


## TODO: repliate the main algorithm here

## (3) process the fitted tree

# (3.1) obtain the (data-defined) support for each of the parameters


# (3.2) obtain the partition --- moment of truth!!


# (3.3) organize all data into single data frame (see partition.R for format)


# ------------------------------------------------------------------------------




## TODO:perform algorithm for a batch approximation (averaged)


## test/validate lambda (closed form vs. numerical) 

## (4) begin main algorithm 
n_partitions = nrow(u_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)








# testing ----------------------------------------------------------------------


u = c(beta, sigmasq)

log_mvnig(u, b_0, V_0, r_0, s_0)

grad(psi_mvn, u, y = y, mu_beta = b_0, V_beta = V_0, a = r_0, b = s_0)

psi_mvn(u, y = y, mu_beta = b_0, V_beta = V_0, a = r_0, b = s_0)




## TODO: check tree output for paramters of the form [beta, sigmasq],
#        where we rely on the order of the first d betas for later parts
#        of the algorithm, e.g., u[1:d], u[d+1]

## thought 1: I think before getting sent into the tree, we already pre-label
## the parameters as u1,...,ud so that we have a handle on every parameter



## TODO: try fitting tree for d' = 3, (beta1, beta2, sigmasq)
## thought 1: previous design just saw me manually creating u_df, i.e,
#  u1 = , u2 = , ... , up = . --- need a way to  avoid doing this, because i 
#  also then compute the closed form integral manually too -> O(d)










# ------------------------------------------------------------------------------














# old stuff --------------------------------------------------------------------


## log of the multivariate normal - inverse gamma density -- 
## log NIG(beta, sigmasq | mu_beta, V_beta, a, b)
log_mvnig = function(beta, sigmasq, mu_beta, V_beta, a, b, d = length(beta)) {
    
    # p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
    #     sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
    #                                                solve(V_beta) %*% 
    #                                                (beta - mu_beta)))
    
    a * log(b) - d / 2 * log(2 * pi) - 0.5 * log_det(V_beta) - lgamma(a) -
        (a + d / 2 + 1) * log(sigmasq) - 
        1 / sigmasq * (b + 0.5 * t(beta - mu_beta) %*% solve(V_beta) %*% 
                           (beta - mu_beta))
}



mvnig = function(beta, sigmasq, mu_beta, V_beta, a, b, d = length(beta)) {
    
    
    p = b^a / ((2 * pi)^(d/2) * sqrt(det(V_beta)) * gamma(a)) * 
        sigmasq^(-a-d/2-1) * exp(-1/sigmasq * (b + 0.5 * t(beta - mu_beta) %*% 
                                                   solve(V_beta) %*% 
                                                   (beta - mu_beta)))
    
    return(p)
}



psi_mvn = function(u, y, n = length(y), d = length(u)) {
    
    loglik = dmvnorm(c(y), mean = X %*% u[1:(d-1)], sigma = u[d] * I_N, log = T)
    logprior = dmvnorm(u[1:(d-1)], mean = b_0, sigma = u[d] * I_D, log = T) + 
        log(dinvgamma(u[d], shape = r_0, scale = s_0))
    
    # print(loglik)
    # print(logprior)
    
    ## TODO: write out log(invgamma) density in closed form to prevent underflow
    
    
    -loglik - logprior
}




