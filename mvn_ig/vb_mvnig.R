# load this LAST to overwrite def preprocess()
# source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 
source("/mnt/N1527194/mlike_approx/mvn_ig/mvn_ig_helper.R")
# rstan_options(auto_write = TRUE)

preprocess = function(u_samps, D, prior) {
    
    #### 2/10 update --- switch over to the psi that we should actually be using 
    # psi_u = apply(u_post, 1, psi_true_mvn, post = post) %>% unname() # (J x 1)
    psi_u = apply(u_samps, 1, psi, prior = prior) %>% unname() # (J x 1)
    
    # (1.2) construct u_df -- this will require some automation for colnames
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(u_samps, psi_u) # (J * N_approx) x (D + 1)
    
    # rename columns 
    names(u_df) = u_df_names
    
    
    return(u_df)
    
}


J         = 10000          # number of MC samples per approximation
N_approx  = 1             # number of approximations
# burn_in   = 2000          # number of burn in draws
# n_chains  = 4             # number of markov chains to run
# stan_seed = 123           # seed

# J_iter = 1 / n_chains * N_approx * J + burn_in 


# K_sims = 1               # num of simulations to run FOR EACH N in N_vec
D = 21
N = 50 # for testing -- comment this line to perform ext. analysis


# set.seed(123)
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

## form the approximation
# post_dat = list(p = p,
#                 a_n = a_n, b_n = b_n, 
#                 mu_star = c(mu_star), V_star = V_star)
# 
# 
# ## start timing here
# start_time <- Sys.time()
# mvnig_fit = stan(file   = 'C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_sampler.stan', 
#                  data    = post_dat,
#                  iter    = J_iter,
#                  warmup  = burn_in,
#                  chains  = n_chains,
#                  refresh = 0) # should give us J * N_approx draws
# 
# # use special preprocess b/c we call psi_true() 
# u_df = preprocess(mvnig_fit, D, post, prior)
# end_time <- Sys.time()
# ## end timing here
# 
# end_time - start_time

### sample from using base-R functions
## start timing here
# start_time <- Sys.time()
sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
beta_samp = matrix(0, J, p)
for (j in 1:J) {
    beta_samp[j,] = rmvnorm(1, mean = mu_star, 
                            sigma =  sigmasq_samp[j] * V_star)
}
u_samps_tmp = data.frame(cbind(beta_samp, sigmasq_samp))
u_df_tmp = preprocess(u_samps_tmp, D, prior)
# end_time <- Sys.time()
## end timing here

# end_time - start_time

hist(u_df_tmp[,2])



# sample from true posterior for sigmasq ~ IG(a_n, b_n)






# sample from approximate posterior for beta ~ N(mu_n, sigmasq_0 * Sigma)
# sigma_0 is the posterior mean of sigmasq
# Sigma is the (scaled) posterior covariance for beta
sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
beta_samp = rmvnorm(J, mean = mu_star, sigma = b_n / (a_n - 1) * V_star)

# beta_samp = matrix(0, nrow = J, ncol = 2)
# for (i in 1:J) {
#     beta_samp[i,] = rmvnorm(1, mean = mu_star, sigma = b_n / (a_n - 1) * V_star)
# }

u_samps = data.frame(cbind(beta_samp, sigmasq_samp))
u_df_mf = preprocess(u_samps, D, prior)
u_df_mf %>% head

data.frame(colMeans(u_df_tmp),
           colMeans(u_df_mf))

hist(u_df_mf[,1])

hml_approx = hml_const(1, D, u_df_tmp, J, prior) 
hml_approx_vb = hml_const(1, D, u_df_mf, J, prior) 

hml_approx$const_vec         # -256.761
hml_approx_vb$const_vec      # -256.761
lil(y, X, prior, post)       # -256.7659


G = 200
approx        = numeric(G)
approx_vb     = numeric(G)
approx_vb_ss  = numeric(G)
true_logml    = numeric(G)

# d_0 = p / 2

# options(warn=0)
for (g in 1:G) {
    
    #### generate data
    X = matrix(rnorm(N * p), N, p) # (N x p) design matrix
    
    eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))
    
    y = X %*% beta + eps # (N x 1) response vector
    
    
    #### compute posterior parameters ------------------------------------
    V_beta_inv = solve(V_beta)
    V_star_inv = t(X) %*% X + V_beta_inv
    
    V_star  = solve(V_star_inv)                                # (p x p)
    mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
    a_n =  a_0 + N / 2 
    b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                              t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                              t(mu_star) %*% V_star_inv %*% mu_star))
    
    # compute MLE estimates for mean, variance of the regression model
    # ybar = X %*% mu_star
    # sigmasq_mle = 1 / N * sum((y - ybar)^2)
    
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
    
    #### TRUE POSTERIOR SAMPLES
    # sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
    # beta_samp = matrix(0, J, p)
    # for (j in 1:J) {
    #     beta_samp[j,] = rmvnorm(1, mean = mu_star, 
    #                             sigma =  sigmasq_samp[j] * V_star)
    #     
    # }
    # u_samps = data.frame(cbind(beta_samp, sigmasq_samp))
    # u_df = preprocess(u_samps, D, prior)
    
    
    #### MFVB SAMPLES
    # sample from approximate posterior for beta ~ N(mu_n, sigmasq_0 * Sigma)
    # sigma_0 is the posterior mean of sigmasq
    # Sigma is the (scaled) posterior covariance for beta
    sigmasq_samp = rinvgamma(J, shape = a_n, scale = b_n)
    beta_samp_1 = rmvnorm(J, mean = mu_star[1:5], 
                          sigma = b_n / (a_n - 1) * V_star[1:5,1:5])
    beta_samp_2 = rmvnorm(J, mean = mu_star[6:10], 
                          sigma = b_n / (a_n - 1) * V_star[6:10,6:10])
    beta_samp_3 = rmvnorm(J, mean = mu_star[11:15], 
                          sigma = b_n / (a_n - 1) * V_star[11:15,11:15])
    beta_samp_4 = rmvnorm(J, mean = mu_star[16:20], 
                          sigma = b_n / (a_n - 1) * V_star[16:20,16:20])
    
    
    beta_samp = cbind(beta_samp_1, beta_samp_2, beta_samp_3, beta_samp_4)
    
    u_samps = data.frame(cbind(beta_samp, sigmasq_samp))
    u_df_mf = preprocess(u_samps, D, prior)
    
    # u_df_mf %>% head
    # data.frame(colMeans(u_df),
    #            colMeans(u_df_mf))
    
    #### compute hybrid approximations for both set of samples
    # hml_approx = hml_const(1, D, u_df, J, prior) 
    hml_approx_vb = hml_const(1, D, u_df_mf, J, prior) 
    
    # approx[g] = hml_approx$const_vec
    approx_vb[g] = hml_approx_vb$const_vec  
    og_part = hml_approx_vb$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    
    ss_fit = fit_resid(og_part, D, 5, prior)
    approx_vb_ss[g] = log_sum_exp(unlist(compute_expterms(ss_fit, D)))
    
    true_logml[g] = lil(y, X, prior, post)      
    
    if (g %% 10 == 0) {
        print(paste("iter: ", g, " -- ", 
                    mean(c(approx_vb[g], approx_vb_ss[g])),
                    ' (', true_logml[g], ')', sep = ''))
    }
}


mean(approx)
mean(approx_vb)
mean(true_logml)

plot(approx, approx_vb)
abline(0, 1)

hist(u_df[,5])