

library(TruncatedNormal)
library(tmg)
library(mvtnorm)


setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")

sourceCpp("C:/Users/ericc/mlike_approx/speedup/trunc_psi.cpp")



set.seed(123)
D = 20
N = 500
I_D = diag(1, D)

n_samps = 10
J       = 5000
B       = 100 # number of replications

hyb_fs  = numeric(B) # store harmonic mean estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator
hyb     = numeric(B) # store harmonic mean estiator


#### generate data -------------------------------------------------------------

# prior mean
mu_0 = rep(0, D)
tau     = 1 / 4          # precision: inverse of variance
sigmasq = 4              # true variance (1 x 1) 

# true value of beta
set.seed(1)
beta = sample(0:10, D, replace = T)
beta = runif(D)
beta = c(runif(D-1, 0, 1), 0)

# generate the regression data -------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector

# compute posterior parameters -------------------------------------------------
Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv =  solve(Q_beta)
b          =  1 / sigmasq * t(X) %*% y
mu_beta    =  Q_beta_inv %*% b


# create prior, post objects to be passed into the hml algorithm 
prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
             Q_beta = Q_beta, b = b)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)



#### compute true log marginal likelihood --------------------------------------

lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
    0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
    1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)

(LIL  = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))

# samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
# u_df = preprocess(data.frame(samples), D, prior)
# hml_approx = hml_const(1, D, u_df, J, prior)
# hml_approx$const_vec


B = 100
hyb = numeric(B)
set.seed(1)

for (i in 1:10) {
    
    # sample from posterior
    samples = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
    u_df = preprocess(data.frame(samples), D, prior)
    
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
        file = 'trunc_d20.RData')
mvnig_d20 = readRDS('trunc_d20.RData')











