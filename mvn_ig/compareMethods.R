

library(rstan)
library(rstudioapi) # running  RStan in parallel via Rstudio
library(reshape2)

rstan_options(auto_write = TRUE)

# DELL_PATH = "C:/Users/chuu/mlike_approx"
LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
setwd(LEN_PATH)

source("partition/partition.R")
source("extractPartition.R")
source("hybrid_approx.R")
source("mvn_ig/mvn_ig_helper.R") # load this LAST to overwrite def preprocess()


# before loading the hme() function, user must have provided a definition for
# likelihood function
source("misc.R") # load the hme() function

N_approx  = 1            # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed
K_sims    = 1            # num of simulations to run FOR EACH N in N_vec

D_vec   = c(3)
J_samps = c(40)
N = 1000 # fix this value for now

set.seed(123)
B = 20

# compute the TRUE log marginal likelihood (fixed for all the B batches) -------

D       = 3                               # dimension of parameter space
p       = D - 1                           # dimension of beta
mu_beta = rep(0, p)                       # prior mean for beta
V_beta  = diag(1, p)                      # scaled precision matrix for beta
a_0     = 2 / 2                           # shape param for sigmasq
b_0     = 1 / 2                           # scale param 
beta    = sample(-10:10, p, replace = T)  # true beta
sigmasq = 4                               # true variance (1 x 1) 
I_p     = diag(1, p)                      # (p x p) identity matrix

## simulate N data points + sample from posterior ------------------------------
X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))

y = X %*% beta + eps # (N x 1) response vector
# ------------------------------------------------------------------------------

## compute posterior parameters ------------------------------------------------
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


J      = 40
J_iter = 1 / n_chains * N_approx * J + burn_in 

## set up the STAN sampler
post_dat = list(p = p,
                a_n = a_n, b_n = b_n, 
                mu_star = c(mu_star), V_star = V_star)

# mvnig_fit = stan(file    = 'mvn_ig/mvn_ig_sampler.stan', 
#                  data    = post_dat,
#                  iter    = J_iter,
#                  warmup  = burn_in,
#                  chains  = n_chains,
#                  refresh = 0) # should give us J * N_approx draws

# true value of the log marginal likelihood:
LIL_mvnig = lil(y, X, prior, post) # -429.3378

# ------------------------------------------------------------------------------



lil_came = numeric(B) # CAME estimator
lil_hme  = numeric(B) # store the log integrated likelihood for each batch
lil_hml  = numeric(B) # store the log integrated likelihood for each batch
# lil_hme_0 = numeric(B)


# number of samples that the corrected arithmetic mean sampler draws from 
# the importance distribution
s_samps = J * B
s_iter = 1 / n_chains * N_approx * s_samps + burn_in 

mvnig_fit = stan(file    = 'mvn_ig/mvn_ig_sampler.stan', 
                 data    = post_dat,
                 iter    = s_iter, ### here
                 warmup  = burn_in,
                 chains  = n_chains,
                 refresh = 0) # should give us J * N_approx draws

u_importance = rstan::extract(mvnig_fit, pars = c("beta", "sigmasq"), 
                              permuted = TRUE)

beta_all    = u_importance$beta %>% data.frame()
sigmasq_all = u_importance$sigmasq 



set.seed(1)
for (b in 1:B) {
    
    print(paste("iter: ", b, "/", B, sep = ""))
    
    # # (0) sample from mu | sigma_sq, y
    # mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    # 
    # # (1) sample from sigma_sq | y
    # sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
    # sample from posterior distribution (via stan)
    
    # NOTE: the following initialization must be declared EVERY TIME otherwise
    # the samples extracted will be exactly the same
    
    # TODO: move this outside later and just sample J * B samples all at once
    # and use the same indexing we use below
    mvnig_fit = stan(file    = 'mvn_ig/mvn_ig_sampler.stan', 
                     data    = post_dat,
                     iter    = J_iter,
                     warmup  = burn_in,
                     chains  = n_chains,
                     refresh = 0) # should give us J * N_approx draws
    
    u_df = preprocess(mvnig_fit, D, post, prior)
    
    
    # (1) compute hybrid app
    hml_approx = hml(N_approx, D, u_df, J, prior) 
    
    lil_hml[b] = hml_approx$hybrid_vec
    lil_hme[b] = hme(u_df, prior, J, D, N)
    
    
    
    # (2) compute corrected arithmetic mean estimator  -------------------------
    
    # this should be able to be done using extractParamSupport
    # A_mu = c(min(mu_post), max(mu_post))
    # A_sigmasq = c(min(sigma_sq_post), max(sigma_sq_post))
     #A_theta = extractSupport(u_df, D) # 1-d intervals stored row-wise
    
    
    ## TODO:  draw from the importance density N-IG
    # mu_s      = rnorm(K, m_n, sqrt(sigma_sq / w_n))
    # sigmasq_s = MCMCpack::rinvgamma(K, shape = r_n / 2, scale = s_n / 2)
    # 
    # mvnig_fit = stan(file    = 'mvn_ig/mvn_ig_sampler.stan', 
    #                  data    = post_dat,
    #                  iter    = s_iter, ### here
    #                  warmup  = burn_in,
    #                  chains  = n_chains,
    #                  refresh = 0) # should give us J * N_approx draws
    # 
    # u_importance = rstan::extract(mvnig_fit, pars = c("beta", "sigmasq"), 
    #                               permuted = TRUE)
    
    start = ((b - 1) * J + 1)
    end = start + J - 1
    beta_imp    = beta_all[start:end,]
    sigmasq_imp = sigmasq_all[start:end]
    
    s_theta = numeric(J)
    lik_j = numeric(J)
    p_theta = numeric(J)
    ind_A = numeric(J)
    
    for (j in 1:J) {
        
        # check containment in A_theta first, if not contained, skip to next
        # iteration of the loop and save on computation
        # 1_A (theta)
        # ind_A = (mu_s >=  A_mu[1] & mu_s <= A_mu[2]) &
        #     (sigmasq_s >= A_sigmasq[1] & sigmasq_s <= A_sigmasq[2]) 
        u_j = c(beta_imp[j,] %>% unlist %>% unname, sigmasq_imp[j])
        ind_A[j] = !sum(!(u_j >= A_theta[,1] & u_j <= A_theta[,2]))
        # if (ind_A == 0) {
        #     next 
        # }
        
        ## TODO: compute 1/s(theta) -- (K x 1) vector of evaluated densities
        # s_theta = dnorm(mu_s, m_n, sqrt(sigma_sq / w_n)) * 
        #     MCMCpack::dinvgamma(sigmasq_s, shape = r_n / 2, scale = s_n / 2)
        
        s_theta[j] = dmvnorm(beta_imp[j,] %>% unlist %>% unname, 
                                mu_star, sigmasq_imp[j] * V_star) * 
            MCMCpack::dinvgamma(sigmasq_imp[j], shape = a_n, scale = b_n)
        
        
        ## TODO: compute prior density
        # p_theta = dnorm(mu_s, m_0, sqrt(sigma_sq / w_0)) * 
        #     MCMCpack::dinvgamma(sigmasq_s, shape = r_0 / 2, scale = s_0 / 2)
        # 
        p_theta[j] = prod(dnorm(beta_imp[j,] %>% unlist %>% unname, 
                             mu_beta, sqrt(sigmasq_imp[j]))) * 
            MCMCpack::dinvgamma(sigmasq_imp[j], shape = a_0, scale = b_0)
        
        ## TODO: compute likelihood
        # lik_q = numeric(J)
        # for (q in 1:J) {
        #     lik_q[q] = prod(dnorm(y, mean = mu_s[q], sd = sqrt(sigmasq_s[q])))
        # }
        #
        lik_j[j] = prod(dnorm(y, mean = X %*% (beta_imp[j,] %>% t), 
                              sqrt(sigmasq_imp[j])))
        
    }
    
    # ind_A
    lil_came[b] = log(mean(1/s_theta * p_theta * lik_j * ind_A))

    
    ## -------------------------------------------------------------------------
    # 
    
} # end of simulation outer loop


hme_df = data.frame(mcmc = 1:B, hml = lil_hml, hme = lil_hme, came = lil_came)

# mean average error (AE, true - estimated)
mean(LIL_mvnig - lil_came)
mean(LIL_mvnig - lil_hml)

# root mean squared error (RMSE)
sqrt(mean((LIL_mvnig - lil_came)^2))
sqrt(mean((LIL_mvnig - lil_hml)^2))

hme_df_long = melt(hme_df, id.vars = "mcmc")

x11()
ggplot(hme_df_long, aes(x = mcmc, y = value, col = variable)) + geom_point() +
    geom_hline(aes(yintercept = LIL_mvnig), linetype = 'dashed', size = 0.9)



# ------------------------------------------------------------------------------



# compute approx over grid of J


# perform analysis for a grid of J
B = 100
J_vec   = c(100, 500, 1000, 2000)  # number of MC samples

lil_hml = data.frame(matrix(0, B, length(J_vec)))
names(lil_hml) = paste('j', J_vec, sep = '')

for (i in 1:length(J_vec)) {
    
    J = J_vec[i]
    J_iter = 1 / n_chains * N_approx * J * B + burn_in 
    
    mvnig_fit = stan(file    = 'mvn_ig/mvn_ig_sampler.stan', 
                     data    = post_dat,
                     iter    = J_iter,
                     warmup  = burn_in,
                     chains  = n_chains,
                     refresh = 0) # should give us J * N_approx draws
    
    print(paste("iter = ", i, ", J = ", J, sep = ''))
    u_df = preprocess(mvnig_fit, D, post, prior)
    
    for (b_i in 1:B) {
        
        # if (b_i %% 20 == 0) print(paste("iter: ", b_i, "/", B, sep = ""))
        
        start = ((b_i - 1) * J + 1)
        end = start + J - 1
        
        u_df_b = u_df[start:end,]
        
        # (1) compute hybrid app
        hml_approx = hml(N_approx, D, u_df_b, J, prior) 
        lil_hml[b_i,i] = hml_approx$hybrid_vec
        
    } # end of simulation outer loop
    
}



## plot results
hme_df = cbind(mcmc = 1:B, lil_hml)

hme_df %>% head
# # mean average error (AE, true - estimated)
# (MAE = round(mean(logZ_mvn - lil_hml), 3))
# 
# # root mean squared error (RMSE)
# (RMSE = round(sqrt(mean((logZ_mvn - lil_hml)^2)), 3))

hme_df_long = melt(hme_df, id.vars = 'mcmc', value.name = 'logZ')


p1 = ggplot(hme_df, aes(x = mcmc, y = j100)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = LIL_mvnig), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(LIL_mvnig, 3), ", D = ", D, 
                       ", J = ", J_vec[1], 
                       ", approx = ", round(mean(hme_df$j100), 2), sep = '')) + 
    ylim(-431, -427) + 
    geom_hline(aes(yintercept = mean(hme_df$j100)), 
               col = 'red', linetype = 'dotdash', size = 1.3)


p2 = ggplot(hme_df, aes(x = mcmc, y = j500)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = LIL_mvnig), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(LIL_mvnig, 3), ", D = ", D, 
                       ", J = ", J_vec[2], 
                       ", approx = ", round(mean(hme_df$j500), 2), sep = '')) + 
    ylim(-431, -427) + 
    geom_hline(aes(yintercept = mean(hme_df$j500)), 
               col = 'red', linetype = 'dotdash', size = 1.3)


p3 = ggplot(hme_df, aes(x = mcmc, y = j1000)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = LIL_mvnig), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(LIL_mvnig, 3), ", D = ", D, 
                       ", J = ", J_vec[3], 
                       ", approx = ", round(mean(hme_df$j1000), 2), sep = '')) + 
    ylim(-431, -427) + 
    geom_hline(aes(yintercept = mean(hme_df$j1000)), 
               col = 'red', linetype = 'dotdash', size = 1.3)

p4 = ggplot(hme_df, aes(x = mcmc, y = j2000)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = LIL_mvnig), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = '',
         title = paste("logZ = ", round(LIL_mvnig, 3), ", D = ", D, 
                       ", J = ", J_vec[4], 
                       ", approx = ", round(mean(hme_df$j2000), 2), sep = '')) + 
    ylim(-431, -427) + 
    geom_hline(aes(yintercept = mean(hme_df$j2000)), 
               col = 'red', linetype = 'dotdash', size = 1.3)

multiplot(p1, p2, p3, p4, cols = 4)






































