
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 


library(Rcpp)
library(RcppEigen)
sourceCpp("C:/Users/ericc/mlike_approx/fast_psi.cpp")

# library(MLmetrics)

# load this LAST to overwrite def preprocess()


J = 10000          # number of MC samples per approximation
D = 20
N = 100 # for testing -- comment this line to perform ext. analysis
n_samps = 10

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


# sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
# beta_post = matrix(0, J, p)
# for (j in 1:J) {
#     beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
# }
sample_beta = function(s2, post) {
    rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
# uncomment line below for 2d
# beta_mat = sapply(sigmasq_post, sample_beta, post = post)

u_samp = data.frame(beta_mat, sigmasq_post)
u_df = preprocess(u_samp, D, prior)

u_df_fast = preprocess(u_samp, D, prior)



# refactored version
hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx = hml_const(1, D, u_df_fast, J, prior) 

hml_approx$param_out %>% 
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs))

hml_approx$const_vec      # -256.761

lil(y, X, prior, post)    # -256.7659


B = 50 # number of replications

hyb_wt1  = numeric(B)  # store harmonic mean estiator
hyb_wt2  = numeric(B)  # store harmonic mean estiator
hyb_wt3  = numeric(B)  # store harmonic mean estiator
# hme     = numeric(B) # store hybrid estimator
# came    = numeric(B) # corrected arithmetic mean estimator
# came_0  = numeric(B)
# bridge  = numeric(B) # bridge estimator (normal)

# other estimators: chib's, bridge, more recent ones?

## bridge sampling specific things ## ------------------------------------------
log_density = function(u, data) {
    -psi(u, data)
}

lb <- c(rep(-Inf, p), 0)
ub <- c(rep(Inf, p), Inf)
colnames(u_samp) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_samp)
## bridge sampling specific things ## ------------------------------------------


# B = 100
# hme     = numeric(B) # store hybrid estimator
# hyb_fs  = numeric(B) # store harmonic mean estiator
# bridge  = numeric(B) # bridge estimator (normal)
set.seed(1)
for (b_i in 1:10) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    u_samp = data.frame(beta_mat, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)
    
    #### (1) hybrid estimator
    # hml_approx = hml_const(1, D, u_df, J, prior)
    # 
    # hml_approx$const_vec
    # source("setup.R") 
    # set.seed(1)
    hybrid = logml(D, u_df, J, prior)
    
    if (any(is.na(hybrid))) {print(paste("error in iteration", b_i)); next;}
    
    # hybrid$all_approx
    # # hybrid$u_df_ss %>% head
    # # hybrid$psi_cand_func %>% head
    # 
    # hybrid$approx_error
    # hybrid$wt_approx1    # error / (2D)
    # hybrid$wt_approx2    # error / sqrt(D)
    # hybrid$wt_approx3
    
    # hybrid$psi_approx %>% head
    # hybrid$psi_true %>% head
    # hybrid$approx_cand
    # 
    # hybrid$approx_cand %>% head
    
    wt_approx = compute_error(hybrid, D)
    wt_approx$all_approx
    wt_approx$wt_approx1
    wt_approx$wt_approx2
    wt_approx$wt_approx3
    
    # 
    # lil(y, X, prior, post)    # -256.7659
    # 
    # hybrid$wt_approx1

    hyb_wt1[b_i] = hybrid$wt_approx1
    hyb_wt2[b_i] = hybrid$wt_approx2
    hyb_wt3[b_i] = hybrid$wt_approx3

    #### (3) corrected arithmetic mean estimator (IS)
    came_result = came_approx(u_df, hml_approx, prior, post, J, D)
    came[b_i] = came_result[1]
    came_0[b_i] = came_result[2]
    
    #### (4) bridge sampling estimator
    # u_samp = as.matrix(u_samp)
    # colnames(u_samp) = names(u_df)[1:D]
    # bridge_result = bridge_sampler(samples = u_samp, 
    #                                log_posterior = log_density,
    #                                data = prior, lb = lb, ub = ub, silent = TRUE)
    # bridge[b_i] = bridge_result$logml
    # if (b_i %% 5 == 0) { 
    #     print(paste("iter ", b_i, ': ',
    #                 "hybrid_wt1 = ", round(mean(hyb_wt1[1:b_i]), 3), '; ',
    #                 "hybrid_wt2 = ", round(mean(hyb_wt2[1:b_i]), 3), '; ',
    #                 "hybrid_avg = ", round(mean(hyb_avg[1:b_i]), 3), '; ',
    #                 sep = '')) 
    # }
    print(paste("iter ", b_i, ': ',
                "hybrid_wt1 = ", round(mean(hyb_wt1[hyb_wt1!=0]), 3), '; ',
                "hybrid_wt2 = ", round(mean(hyb_wt2[hyb_wt2!=0]), 3), '; ',
                "hybrid_wt3 = ", round(mean(hyb_wt3[hyb_wt3!=0]), 3), '; ',
                sep = '')) 
}

LIL = lil(y, X, prior, post)
approx = data.frame(LIL, 
                    hyb_wt1 = hyb_wt1[hyb_wt1!=0], 
                    hyb_wt2 = hyb_wt2[hyb_wt2!=0], 
                    hyb_wt3 = hyb_wt3[hyb_wt3!=0])
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)








