

###

source("partition.R")

library(dplyr)

set.seed(123)

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N = 50

# generate 50 samples from N(mu, sigma_sq)
y = rnorm(N, mu, sqrt(sigma_sq))

ybar = mean(y)

# compute posterior parameters
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2

p_y = (pi)^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
    s_0^(r_0 / 2) / s_n^(r_n / 2)


(LIL = log(p_y)) # -113.143 (paper says -117, but difference arises from RNG)


# re-eval by first modifying rpart() parameters to encourage a deeper tree

# normal inverse gamma density with params (mu, lambda, alpha, beta)
nig = function(x, sigmasq, mu, lambda, alpha, beta) {
    
    sqrt(lambda) / (sqrt(sigmasq * 2 * pi)) * beta^alpha / gamma(alpha) * 
        sigmasq^(-alpha - 1) * 
        exp(-(2 * beta + lambda * (x - mu)^2) / (2 * sigmasq))
    
}

psi_true = function(mu, sigma_sq, m_n, w_n, r_n, s_n) {
    
    # p (mu, sigma_sq | y ) 
    #    = N (mu | m_n, sigma_sq / w_n) * IG (sigma_sq | r_n / 2, s_n / 2)
    #    = NIG (mu, sigma_sq | )
    
    log_mu_pdf   = dnorm(mu, m_n, sqrt(sigma_sq / w_n), log = T)
    log_sigma_sq = log(MCMCpack::dinvgamma(sigma_sq, r_n / 2, s_n / 2))
    log_p_mu_sigmasq = log_mu_pdf + log_sigma_sq
    
    out = -log_p_mu_sigmasq

    return(-log(nig(mu, sigma_sq, m_n, w_n, r_n / 2, s_n / 2)))
}

# psi_tilde(u) = c_k + lambda_k'u
psi = function(mu, sigmasq, y, m_0, w_0, r_0, s_0) {
    loglik = sum(dnorm(y, mu, sqrt(sigmasq), log = TRUE))
    logprior = log(nig(mu, sigmasq, m_0, w_0, r_0 / 2, s_0 / 2))
    
    out = -loglik - logprior
    # out = -logprior
    
    return(out)
}

RECURSIVE = TRUE
N_iters = 100
def_approx = numeric(N_iters)
rec_approx = numeric(N_iters)
for (t in 1:N_iters) {
    print(paste("iter", t))
    # generate samples from the posterior probability to form the HME estimator
    J = 5000 # number of random draws used per estimate
    
    # (0) sample from mu | sigma_sq, y
    mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
    
    # (1) sample from sigma_sq | y
    sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
    
    psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)
    
    # input for paramPartition() MUST have parameter names u1, u2, ... up
    u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u)    # (J x 3)
    
    # fit decision tree
    ### use rpart to fit partition
    nig_rpart = rpart(psi_u ~ ., u_df)
    
    # plot the tree (matches general structure of the tree returned from tree())
    # plot(nig_rpart)
    # text(nig_rpart, cex = 0.7)
    
    ### obtain partition
    nig_support = rbind(c(min(mu_post), max(mu_post)),
                        c(min(sigma_sq_post), max(sigma_sq_post)))
    nig_partition = paramPartition(nig_rpart, nig_support)  # partition.R
    
    # organize all data into single data frame --> ready for approximation
    param_out = u_star(nig_rpart, u_df, nig_partition)
    
    if (RECURSIVE) {
        n_partitions = nrow(nig_partition)
        c_k = numeric(n_partitions)
        zhat = numeric(n_partitions)
        for (k in 1:n_partitions) {
            
            # u_star_k = (mu_k, sigma_sq_k)
            c_k[k] = exp(-psi_hat(param_out[k,]$u1_star, 
                                  param_out[k,]$u2_star,
                                  y, m_0, w_0, r_0, s_0)) # (1 x 1)
            
            l_k = lambda(param_out[k,]$u1_star, param_out[k,]$u2_star,
                         y, m_0, w_0, r_0, s_0)
            
            # 1st param calculation
            p1 = -1/l_k[1] * # exp(l_k[1] * param_out[k,]$u1_star) *
                exp(-l_k[1] * (param_out[k,]$u1_ub - param_out[k,]$u1_lb))
            
            # 2nd param calculation
            p2 = -1/l_k[2] * # exp(l_k[2] * param_out[k,]$u2_star) *
                exp(-l_k[2] * (param_out[k,]$u2_ub - param_out[k,]$u2_lb))
            
            zhat[k] = c_k[k] * p1 * p2
        }
    }
    
    def_approx[t] = log(sum(zhat))
        
    u_df_leaf = u_df %>% mutate(leaf_id = nig_rpart$where)  
    leaf_id_max = as.numeric(names(sort(table(u_df_leaf$leaf_id), 
                                        decreasing = T)))[1]
    
    
    u_df_sub = u_df_leaf %>% filter(leaf_id == leaf_id_max)
    
    rpart_sub = rpart(psi_u ~ u1 + u2, u_df_sub)
    
    sub_support = with(u_df_sub, rbind(c(min(u1), max(u1)),
                                       c(min(u2), max(u2))))
    
    sub_partition = paramPartition(rpart_sub, sub_support)  # partition.R
    
    # organize all data into single data frame --> ready for approximation
    param_out_sub = u_star(rpart_sub, u_df_sub, sub_partition)
    
    param_mod = do.call(rbind, list(param_out %>% filter(leaf_id != leaf_id_max),
                                    param_out_sub))
    
    n_partitions = nrow(param_mod)
    c_k = numeric(n_partitions)
    zhat_mod = numeric(n_partitions)
    for (k in 1:n_partitions) {
        
        # print(paste("iter", k))
        # u_star_k = (mu_k, sigma_sq_k)
        c_k[k] = exp(-psi_hat(param_mod[k,]$u1_star, 
                              param_mod[k,]$u2_star,
                              y, m_0, w_0, r_0, s_0)) # (1 x 1)
        
        l_k = lambda(param_mod[k,]$u1_star, param_mod[k,]$u2_star,
                     y, m_0, w_0, r_0, s_0)
        
        # 1st param calculation
        p1 = -1/l_k[1] * exp(-l_k[1] * (param_mod[k,]$u1_ub - param_mod[k,]$u1_lb))
        
        # 2nd param calculation
        p2 = -1/l_k[2] * exp(-l_k[2] * (param_mod[k,]$u2_ub - param_mod[k,]$u2_lb))
        
        zhat_mod[k] = c_k[k] * p1 * p2
    }
    
    rec_approx[t] = log(sum(zhat_mod))
}



# after default
approx_default = def_approx
approx_df = data.frame(def = approx_default, lil = LIL) %>% 
    mutate(iter = 1:length(approx_default))
approx_df_full = na.omit(approx_df)
approx_df_full = approx_df_full %>% mutate(iter = 1:nrow(approx_df_full))

# after one round of recursive partitioning
approx_rec = rec_approx
approx_rec_df = data.frame(rec = approx_rec, lil = LIL) %>% 
    mutate(iter = 1:length(approx_rec))

# HME estimator

## run the code in hme_estimator.R

ggplot(approx_df, aes(iter, def)) + geom_point() + 
    geom_point(aes(hme_df$mcmc, hme_df$hme), col = 'green', shape = 5) + 
    geom_point(aes(approx_rec_df$iter, approx_rec_df$rec), col = 'purple', shape = 8) + 
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 1.5, color = "red")


ggplot(approx_df, aes(iter, def)) + geom_point() + 
    geom_point(aes(approx_rec_df$iter, approx_rec_df$rec), col = 'purple', shape = 8) + 
    geom_hline(aes(yintercept = LIL), linetype = 'dashed', size = 1.5, color = "red")


mean(approx_default, na.rm = T)
var(approx_default, na.rm = T)
mean(approx_rec, na.rm = T)
var(approx_rec, na.rm = T)



### use rpart to fit partition
nig_rpart = rpart(psi_u ~ ., u_df, cp = 0.002)



### obtain partition
nig_support = rbind(c(min(mu_post), max(mu_post)),
                    c(min(sigma_sq_post), max(sigma_sq_post)))
nig_partition = paramPartition(nig_rpart, nig_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(nig_rpart, u_df, nig_partition)

n_partitions = nrow(nig_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)
for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    c_k[k] = exp(-psi_hat(param_out[k,]$u1_star, 
                          param_out[k,]$u2_star,
                          y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    l_k = lambda(param_out[k,]$u1_star, param_out[k,]$u2_star,
                 y, m_0, w_0, r_0, s_0)
    
    # 1st param calculation
    p1 = -1/l_k[1] * exp(-l_k[1] * (param_out[k,]$u1_ub - param_out[k,]$u1_lb))
    
    # 2nd param calculation
    p2 = -1/l_k[2] * exp(-l_k[2] * (param_out[k,]$u2_ub - param_out[k,]$u2_lb))
    
    zhat[k] = c_k[k] * p1 * p2
}

log(sum(zhat))



# plot the tree (matches general structure of the tree returned from tree())
plot(nig_rpart)
text(nig_rpart, cex = 0.7)



