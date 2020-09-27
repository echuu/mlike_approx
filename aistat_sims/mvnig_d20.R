
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 


library(Rcpp)
library(RcppEigen)
sourceCpp("C:/Users/ericc/mlike_approx/fast_psi.cpp")

# library(MLmetrics)

# load this LAST to overwrite def preprocess()


J = 10000          # number of MC samples per approximation
D = 25
N = 500


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
prior = list(V_beta = V_beta, mu_beta = mu_beta,  a_0 = a_0, b_0 = b_0,
             y = y, X = X, V_beta_inv = V_beta_inv)
# store posterior parameters
post  = list(V_star =  V_star, mu_star = mu_star, a_n = a_n, b_n = b_n, 
             V_star_inv = V_star_inv)

sample_beta = function(s2, post) { 
    rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

(LIL = lil(y, X, prior, post))    # -256.7659

B = 100
hyb = numeric(B)
set.seed(1)

for (i in 1:B) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    u_samp = data.frame(beta_mat, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)
    
    #### (1) hybrid estimator
    hybrid = hybrid_ml(D, u_df, J, prior)
    
    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}
    
    hyb[i] = hybrid$zhat
    
    avg_hyb = mean(hyb[hyb!=0])
    print(paste("iter ", i, ': ',
                "hybrid = ", round(avg_hyb, 3), '; ',
                "ae = ", LIL - avg_hyb,
                sep = '')) 
}

approx = data.frame(LIL, hyb = hyb[hyb!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)

saveRDS(list(J = J, D = D, N = N, approx_df = approx), 
        file = 'mvnig_d20.RData')
mvnig_d20 = readRDS('mvnig_d20.RData')






