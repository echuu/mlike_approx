



# setup global environment, load in algo functions
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

# load this LAST to overwrite def preprocess()
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
options(mc.cores = parallel::detectCores()) 


J         = 1000          # number of MC samples per approximation
N_approx  = 1             # number of approximations
burn_in   = 2000          # number of burn in draws
n_chains  = 4             # number of markov chains to run
stan_seed = 123           # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 


K_sims = 1               # num of simulations to run FOR EACH N in N_vec
D = 20
N = 200 # for testing -- comment this line to perform ext. analysis


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

n_samps = 10

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

## form the approximation
# post_dat = list(p = p,
#                 a_n = a_n, b_n = b_n,
#                 mu_star = c(mu_star), V_star = V_star)
# 
# mvnig_fit = stan(file   = 'C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_sampler.stan',
#                  data    = post_dat,
#                  iter    = J_iter,
#                  warmup  = burn_in,
#                  chains  = n_chains,
#                  refresh = 0) # should give us J * N_approx draws
# 
# stan_samp = data.frame(rstan::extract(mvnig_fit, pars = c("beta", "sigmasq"),
#                            permuted = TRUE))

# use special preprocess b/c we call psi_true() 

sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_post = matrix(0, J, p)
for (j in 1:J) {
    beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
}
u_samp = data.frame(beta_post, sigmasq_post)
u_df = preprocess(u_samp, D, prior)


# refactored version
hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx$param_out %>% 
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs))

hml_approx$const_vec      # -256.761
# came_approx(u_df, hml_approx, prior, post, J, D)

lil(y, X, prior, post)    # -256.7659

abs(hml_approx$const_vec - lil(y, X, prior, post))


#### begin simulations ---------------------------------------------------------

B = 1000 # number of replications
J = 5000 # number of MCMC samples per replication

hyb_fs  = numeric(B) # store harmonic mean estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator
hyb     = numeric(B) # store harmonic mean estiator
hme     = numeric(B) # store hybrid estimator
came    = numeric(B) # corrected arithmetic mean estimator

# other estimators: chib's, bridge, more recent ones?

for (b in 1:B) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_post = matrix(0, J, p)
    for (j in 1:J) {
        beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
    }
    u_samp = data.frame(beta_post, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)
    
    
    #### (1) hybrid estimator
    hml_approx = hml_const(1, D, u_df, J, prior)
    og_part = hml_approx$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    ss_part = fit_resid(og_part, D, n_samps, prior)
    ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
    
    hyb_fs[b] = hml_approx$const_vec
    hyb_ss[b] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    hyb_ts[b] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    hyb[b] = mean(c(hyb_fs[b], hyb_ss[b], hyb_ts[b]))
    
    #### (2) harmonic mean estimator
    hme[b] = hme_approx(u_df, prior, J, D, N)
    
    #### (3) corrected arithmetic mean estimator (IS)
    came[b] = came_approx(u_df, hml_approx, prior, post, J, D)
    
    if (b %% 100 == 0) { print(paste("iter:", b)) }
}

LIL = lil(y, X, prior, post)
LIL

x1 = round(mean(hyb_fs), 3)
x2 = round(mean(hyb_ss), 3)
x3 = round(mean(hyb_ts), 3)
mean(hyb)

(mean(c(x1,x2,x3)) + median(c(x1, x2, x3)))/2

pred_df = do.call(cbind, list(hyb_fs, hyb_ss, hyb_ts))
pred_df %>% head

pred_bar = apply(pred_df, 1, mean)
pred_med = apply(pred_df, 1, median)

round(mean((LIL - hyb_fs)), 3)
mean(hyb_fs)



pred_vec = (pred_bar + pred_med)/2
mean(pred_vec)
sd(pred_vec)
round(mean(LIL - pred_vec), 3)
round(sqrt(mean((LIL - pred_vec)^2)), 3)

round(mean(hme), 3)
# mean(ame)
round(sd(came), 3)
round(mean(hyb_fs), 3)
round(mean(hyb_ss), 3)
round(mean(hyb_ts), 3)


mean(hyb_fs)
sd(hyb_fs)

round(mean(LIL - hyb), 3)
round(mean(LIL - hyb_fs), 3)
round(sqrt(mean((LIL - hyb_fs)^2)), 3)

# mean(LIL - ame)
round(mean(LIL - came), 3)

round(sqrt(mean((LIL - hyb)^2)), 3)
round(sqrt(mean((LIL - hme)^2)), 3)
# sqrt(mean((LIL - ame)^2))
round(sqrt(mean((LIL - came)^2)), 3)



# ------------------------------------------------------------------------------

# contains: leaf_id, u1_lb, u1_ub, ... , uD_lb, uD_ub, n_obs
part_0 = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

(orig_partition = hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs)))


K = length(part_set)

# initialize a list to store the vector containing the terms in exponential
# for each of the sub-partitions
# kth elmt is an s_k dim vector of terms that are to be exponentiated
# at the very end, all entries are unlisted and evaluated via log-sum-exp
exp_terms = vector("list", K) 
ck_star_list = vector("list", K)

perc_thresh = sort(orig_partition$perc, decreasing = T)

for (k in 1:K) {
    
    
    PERC_K = orig_partition[k,]$perc
    
    if (PERC_K >= perc_thresh[9]) {
        print("using original partition")
        # exp_terms[[k]] = hml_approx$const_approx[k]
        
        N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        part_k = part_0 %>%           # set of lower/upper bounds
            dplyr::filter(leaf_id == part_set[k]) %>%
            dplyr::select(-c(leaf_id, n_obs))

        # sample uniformly from each lower/upper bound pair to form a D-dim vector
        part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)

        resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1],
                                upper = part_k_long[,2]) %>% data.frame

        u_df_k = preprocess_resample(resamp_k, D, prior) # N_k_p x (D_u + 1)

        c_k_approx = hml_const_mod(1, D, u_df_k, N_k_p, prior)

        ck_star_list[[k]] = c_k_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)

        exp_terms[[k]] = c_k_approx$const_approx
        
    } else {
        N_k_p = part_0$n_obs[k] * 10  # number of (re)samples to draw from part k
        part_k = part_0 %>%           # set of lower/upper bounds
            dplyr::filter(leaf_id == part_set[k]) %>% 
            dplyr::select(-c(leaf_id, n_obs))
        
        # sample uniformly from each lower/upper bound pair to form a D-dim vector
        part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
        
        resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1], 
                                upper = part_k_long[,2]) %>% data.frame
        
        u_df_k = preprocess_resample(resamp_k, D, prior) # N_k_p x (D_u + 1)
        
        c_k_approx = hml_const(1, D, u_df_k, N_k_p, prior)
        
        ck_star_list[[k]] = c_k_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) 
        
        exp_terms[[k]] = c_k_approx$const_approx
    }
    
}

all_terms = exp_terms %>% unlist

log_sum_exp(all_terms)    # 
hml_approx$const_vec      # -256.761
lil(y, X, prior, post)    # -256.7659

abs(hml_approx$const_vec - lil(y, X, prior, post))
abs(log_sum_exp(all_terms) - lil(y, X, prior, post))





