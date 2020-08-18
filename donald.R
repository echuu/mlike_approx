

#### example run through


#### LOAD LIBRARIES
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

# load this LAST to overwrite def preprocess()
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 

options(scipen = 999)
options(dplyr.summarise.inform = FALSE)


J         = 1000          # number of MC samples per approximation
N_approx  = 1             # number of approximations
D = 2   # dimension of parameter
N = 200 # sample size

#### DATA GENERATION 
set.seed(123)
p       = D - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
a_0     = 2 / 2            # shape param for sigmasq
b_0     = 1 / 2            # scale param 
beta    = sample(-10:10, p, replace = T)
sigmasq = 4                # true variance (1 x 1) 

I_p = diag(1, p)           # (p x p) identity matrix
I_N = diag(1, N)           # (N x N) identity matrix

X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))

y = X %*% beta + eps # (N x 1) response vector
# ------------------------------------------------------------------

## compute posterior parameters ------------------------------------
V_beta_inv = solve(V_beta)
V_star_inv = t(X) %*% X + V_beta_inv

V_star  = solve(V_star_inv)                                # (p x p)
mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
a_n =  a_0 + N / 2 
b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                          t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                          t(mu_star) %*% V_star_inv %*% mu_star))

# compute MLE estimates for mean, variance of the regression model
ybar = X %*% mu_star
sigmasq_mle = 1 / N * sum((y - ybar)^2)

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


sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_post = matrix(0, J, p)
for (j in 1:J) {
    beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
}
u_samp = data.frame(beta_post, sigmasq_post)


u_df = preprocess(u_samp, D, prior) # evaluate each u_i w/ Psi()



hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx$const_vec


logml(D, u_df, J, prior)














