

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 

library(bridgesampling)
library(Rcpp)
library(RcppEigen)
sourceCpp("C:/Users/ericc/mlike_approx/fast_psi.cpp")
# 
# preprocess = function(u_samps, D, prior) {
#     
#     #### 2/10 update --- switch over to the psi that we should actually be using 
#     # psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)
#     psi_u = apply(u_samps, 1, psi, prior = prior) %>% unname() # (J x 1)
#     
#     # (1.2) construct u_df -- this will require some automation for colnames
#     u_df_names = character(D + 1)
#     for (d in 1:D) {
#         u_df_names[d] = paste("u", d, sep = '')
#     }
#     u_df_names[D + 1] = "psi_u"
#     
#     # populate u_df
#     u_df = cbind(u_samps, psi_u) # (J * N_approx) x (D + 1)
#     
#     # rename columns 
#     names(u_df) = u_df_names
#     
#     return(u_df)
#     
# }


J  = 200          # number of MC samples per approximation

# J_iter = 1 / n_chains * N_approx * J + burn_in 


# K_sims = 1               # num of simulations to run FOR EACH N in N_vec
D = 10
N = 50 # for testing -- comment this line to perform ext. analysis


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

(LIL = lil(y, X, prior, post))    # -256.7659
# sample from true posterior for sigmasq ~ IG(a_n, b_n)
# sample from approximate posterior for beta ~ N(mu_n, sigmasq_0 * Sigma)
# sigma_0 is the posterior mean of sigmasq
# Sigma is the (scaled) posterior covariance for beta
# beta_samp = rmvnorm(J, mean = mu_star, sigma = b_n / (a_n - 1) * V_star)
# 
sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
beta_samp = matrix(0, nrow = J, ncol = 2)
for (i in 1:J) {
    beta_samp[i,] = rmvnorm(1, mean = mu_star, sigma = b_n / (a_n - 1) * V_star)
}
# 
# u_samps = data.frame(cbind(beta_samp, sigmasq_samp))
# u_df_mf = preprocess(u_samps, D, prior)
# u_df_mf %>% head
# 
# data.frame(colMeans(u_df_tmp),
#            colMeans(u_df_mf))
# 
# hist(u_df_mf[,1])
# 
# hml_approx = hml_const(1, D, u_df_tmp, J, prior) 
# hml_approx_vb = hml_const(1, D, u_df_mf, J, prior) 
# 
# hml_approx$const_vec         # -256.761
# hml_approx_vb$const_vec      # -256.761




# J = 100 for MF simulation -> 0.55
# J = 1000 for correct -> 0.65 (want this to be lower)

# J = 2000

log_density = function(u, data) {
    -psi(u, data)
}

J = 1000
B = 100
hyb        = numeric(B)
hyb_mf     = numeric(B)
bridge     = numeric(B)
# d_0 = p / 2

# options(warn=0)
set.seed(1)
for (i in 1:B) {
    
    # #### TRUE POSTERIOR SAMPLES
    # sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    # beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    
    sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
    beta_samp = matrix(0, nrow = J, ncol = p)
    for (j in 1:J) {
        beta_samp[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_samp[j] * V_star)
    }
    u_samps = data.frame(beta_samp, sigmasq_samp)
    u_df = preprocess(u_samps, D, prior)
    hybrid = hybrid_ml(D, u_df, J, prior)
    hybrid$zhat
    hyb[i] = hybrid$zhat
    
    # #### MFVB SAMPLES
    # sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
    # beta_samp_1 = rmvnorm(J, mean = mu_star[1:3],
    #                       sigma = mean(sigmasq_samp) * V_star[1:3,1:3])
    # beta_samp_2 = rmvnorm(J, mean = mu_star[4:6],
    #                       sigma = mean(sigmasq_samp) * V_star[4:6,4:6])
    # beta_samp_3 = rmvnorm(J, mean = mu_star[7:9],
    #                       sigma = mean(sigmasq_samp) * V_star[7:9,7:9])
    # 
    # # beta_samp = cbind(beta_samp_1, beta_samp_2)
    # beta_samp = cbind(beta_samp_1, beta_samp_2, beta_samp_3)
    # 
    # u_samps = data.frame(cbind(beta_samp, sigmasq_samp))
    # u_df = preprocess(u_samps, D, prior)
    # 
    # hybrid_mf = hybrid_ml(D, u_df, J, prior)
    # hybrid_mf$zhat
    # hyb[i] = hybrid_mf$zhat
    
    
    lb <- c(rep(-Inf, p), 0)
    ub <- c(rep(Inf, p), Inf)
    u_samps = as.matrix(u_samps)
    colnames(u_samps) = names(u_df)[1:D]
    names(lb) <- names(ub) <- colnames(u_samps)
    bridge_result = bridge_sampler(samples = u_samps,
                                   log_posterior = log_density, data = prior, 
                                   lb = lb, ub = ub, silent = TRUE)
    bridge_result$logml
    bridge[i] = bridge_result$logml
    
    print(paste("iter ", i, ': ',
                "hyb_mf = ", round(mean(hyb[hyb!=0]), 3),
                '; ', "mae = ", round(mean(abs(LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(mean(bridge[bridge!=0]), 3), '; ',
                "mae = ", round(mean(abs(LIL - bridge[bridge!=0])), 4),
                # "came = ", round(mean(came[came!=0]), 3), '; ', 
                # "mae = ", round(mean(abs(LIL - came[came!=0])), 4),
                sep = '')) 
}


approx = data.frame(LIL, hyb = hyb[hyb!=0], bridge = bridge[bridge!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(abs(LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)

data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans((LIL - approx)),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)



abline(0, 1)

hist(u_df[,5])