setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

# load this LAST to overwrite def preprocess()
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 




J = 1000          # number of MC samples per approximation
D = 60
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

sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
u_samp = data.frame(beta_mat, sigmasq_post)
u_df = preprocess(u_samp, D, prior)

# refactored version
hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx$param_out %>% 
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs))

hml_approx$const_vec      # -256.761
lil(y, X, prior, post)    # -256.7659



B = 100 # number of replications
# J = 5000 # number of MCMC samples per replication

hyb_fs  = numeric(B) # store harmonic mean estiator
hyb_ss  = numeric(B) # store harmonic mean estiator
hyb_ts  = numeric(B) # store harmonic mean estiator
hyb     = numeric(B) # store harmonic mean estiator
hme     = numeric(B) # store hybrid estimator
came    = numeric(B) # corrected arithmetic mean estimator
came_0  = numeric(B)
bridge  = numeric(B) # bridge estimator (normal)

# other estimators: chib's, bridge, more recent ones?

sample_beta = function(s2, post) {
   rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

## bridge sampling specific things ## ------------------------------------------
log_density = function(u, data) {
    -psi(u, data)
}

lb <- c(rep(-Inf, p), 0)
ub <- c(rep(Inf, p), Inf)
colnames(u_samp) = names(u_df)[1:D]
names(lb) <- names(ub) <- colnames(u_samp)
## bridge sampling specific things ## ------------------------------------------


B = 100
hme     = numeric(B) # store hybrid estimator
hyb_fs  = numeric(B) # store harmonic mean estiator
bridge  = numeric(B) # bridge estimator (normal)
set.seed(1)
for (b_i in 1:B) {
    
    # sample from posterior
    sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
    beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
    u_samp = data.frame(beta_mat, sigmasq_post)
    u_df = preprocess(u_samp, D, prior)
    
    #### (1) hybrid estimator
    hml_approx = hml_const(1, D, u_df, J, prior)
    # og_part = hml_approx$param_out %>%
    #     dplyr::select(-c(psi_choice, logQ_cstar))
    # ss_part = fit_resid(og_part, D, n_samps, prior)
    # ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
    # 
    # hyb_fs[b_i] = hml_approx$const_vec
    # hyb_ss[b_i] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    # hyb_ts[b_i] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    # hyb[b_i] = mean(c(hyb_fs[b_i], hyb_ss[b_i], hyb_ts[b_i]))
    
    #### (2) harmonic mean estimator
    # hme[b_i] = hme_approx(u_df, prior, J, D, N)
    
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
    
    if (b_i %% 10 == 0) { 
        print(paste("iter: ", b_i, 
                    # "hyb = ", round(mean(hyb[1:b_i]), 3),
                    "came = ", round(mean(came[1:b_i]), 3),
                    "came0 = ", round(mean(came_0[1:b_i]), 3))) 
    }
}

LIL = lil(y, X, prior, post)    # -256.7659
approx = data.frame(bridge, hme)

approx = data.frame(hyb[1:B], came[1:B], came_0)
data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
           ae = colMeans(LIL - approx),
           rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3)

mean(hyb_fs[1:B])













