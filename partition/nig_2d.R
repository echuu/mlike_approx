

source("partition.R")

library(dplyr)


#### specify hyperparameters
set.seed(1)

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N = 50


#### generate data
y = rnorm(N, mu, sqrt(sigma_sq))

#### compute posterior parameters
ybar = mean(y)
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2


#### compute true (log) marginal likelihood
# p_y = (pi)^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
#   s_0^(r_0 / 2) / s_n^(r_n / 2)

# (LIL = log(p_y)) # -108.877

LIL = log_p_y = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
           lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
           r_n / 2 * log(s_n)

print(LIL)

#### specify prior distribution density
nig = function(x, sigmasq, mu, lambda, alpha, beta) {
  
  sqrt(lambda) / (sqrt(sigmasq * 2 * pi)) * beta^alpha / gamma(alpha) * 
    sigmasq^(-alpha - 1) * 
    exp(-(2 * beta + lambda * (x - mu)^2) / (2 * sigmasq))
  
}


# ------------------------------------------------------------------------------

## psi_true:     the true negative log posterior, not available in practice, but 
##               we can evaluate it
## psi:          the negative log posterior as described in the notes, 
##               = -loglik - logprior 
## psi_tilde:    approximation of psi as described in the notes
##               = c_k + lambda_k'u

## note: (1) in harder problems, we only have psi(), psi_tilde()
##       (2) in conjugate examples, we know the posterior in closed form, so we
##           fit the tree using psi_true(), but the approximation part still 
##           uses psi()

#### (*) specify psi_true (in more difficult problems, this one shouldn't be 
#### avaiable, but we're cheating a little to use this when fitting the 
#### decision tree)

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


#### (**) specify psi = -loglik - logprior
psi = function(mu, sigmasq, y, m_0, w_0, r_0, s_0) {
    loglik = sum(dnorm(y, mu, sqrt(sigmasq), log = TRUE))
    logprior = log(nig(mu, sigmasq, m_0, w_0, r_0 / 2, s_0 / 2))
    
    out = -loglik - logprior

    return(out)
}


#### specify lambda = grad(psi)
lambda = function(mu_star, sigmasq_star, y, m_0, w_0, r_0, s_0) {
    
    n = length(y)
    
    lambda1 = -1 / sigmasq_star * sum((y - mu_star)) + 
        w_0 / sigmasq_star * (mu_star - m_0)
    
    lambda2 = n / (2 * sigmasq_star)  - 
        1 / (2 * sigmasq_star^2) * sum((y - mu_star)^2) +
        (r_0 / 2 + 3 / 2) / sigmasq_star - 
        1 / (2 * sigmasq_star^2) * (w_0 * (mu_star - m_0)^2 + s_0)
    
    return(c(lambda1, lambda2))
    
}

#### (***) specify psi_tilde = c_k + lambda_k'u
#### this is calculated within the main loop -- move it out later


#### specify algorithm settings
N_iters = 100                     # number of approximations
# rec_approx = numeric(N_iters)   # storage for r.p. approximations

approx_lil = function(N_iters, y, y_bar,
                      m_0, w_0, r_0, s_0,
                      m_n, w_n, r_n, s_n) {
    
    #### algorithm: main loop
    
    def_approx = numeric(N_iters)   # storage for default approximations (no r.p.)
    
    for (t in 1:N_iters) {
        
        if (t %% 10 == 0) {
            print(paste("iter", t))
        }
        
        # generate samples from the posterior probability to form the HME estimator
        J = 3000 # number of random draws used per estimate
        
        # (0) sample from mu | sigma_sq, y
        mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)
        
        # (1) sample from sigma_sq | y
        sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)
        
        psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)
        
        # input for paramPartition() MUST have parameter names u1, u2, ... up
        u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u) # (J x 3)
        
        # fit decision tree
        ### use rpart to fit partition
        u_rpart = rpart(psi_u ~ ., u_df)
        
        # plot the tree (matches general structure of the tree returned from tree())
        # plot(nig_rpart)
        # text(nig_rpart, cex = 0.7)
        
        ### obtain partition
        u_support = rbind(c(min(mu_post), max(mu_post)),
                          c(min(sigma_sq_post), max(sigma_sq_post)))
        
        u_partition = paramPartition(u_rpart, u_support)  # partition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition)
        
        n_partitions = nrow(u_partition)
        c_k = numeric(n_partitions)
        zhat = numeric(n_partitions)
        for (k in 1:n_partitions) {
            
            # u_star_k = (mu_k, sigma_sq_k)
            c_k[k] = exp(-psi(param_out[k,]$u1_star, 
                                  param_out[k,]$u2_star,
                                  y, m_0, w_0, r_0, s_0)) # (1 x 1)
            
            l_k = lambda(param_out[k,]$u1_star, param_out[k,]$u2_star,
                         y, m_0, w_0, r_0, s_0)
            
            # 1st param calculation
            p1 = -1 / l_k[1] * 
                exp(-l_k[1] * (param_out[k,]$u1_ub - param_out[k,]$u1_lb))
            
            # 2nd param calculation
            p2 = -1 / l_k[2] * 
                exp(-l_k[2] * (param_out[k,]$u2_ub - param_out[k,]$u2_lb))
            
            zhat[k] = c_k[k] * p1 * p2
        }
        
        def_approx[t] = log(sum(zhat))
    }
    
    return(def_approx)
}


set.seed(1)
def_approx = approx_lil(N_iters, y, y_bar,
                        m_0, w_0, r_0, s_0,
                        m_n, w_n, r_n, s_n)

# evaluate the mean, variance of the approximations
mean(def_approx, na.rm = TRUE) # -106.1804
var(def_approx, na.rm = TRUE)  # 4.491305

# compare approximation with the true value of the log marginal likelihood
approx_default = def_approx
approx_df = data.frame(def = approx_default, lil = LIL) %>% 
    mutate(iter = 1:length(approx_default))

ggplot(approx_df, aes(iter, def)) + geom_point() + 
    geom_hline(aes(yintercept = LIL), 
               linetype = 'dashed', size = 1.5, color = "red")



#### ---------------------------------------------------------------------------

#### log marginal likelihood analysis ------------------------------------------

N_vec = c(50, 100, 200, 500)
N_vec = c(50, 100, 200, 500, 1000, 2000, 10000)

#### generate data
N = 50
y = rnorm(N, mu, sqrt(sigma_sq))

#### compute posterior parameters
ybar = mean(y)
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2


#### compute true (log) marginal likelihood
# p_y underflows for relatively small N -- use log_p_y definition instead
# p_y = (pi)^(-N / 2) * (w_0 / w_n)^(1/2) * gamma(r_n / 2) / gamma(r_0 / 2) * 
#     s_0^(r_0 / 2) / s_n^(r_n / 2)
# (LIL = log(p_y)) # -108.877

(log_p_y = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
    lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
    r_n / 2 * log(s_n))


## plot log-marginal-likelihood vs log n
## overlay approximate log-marginal likelihood vs log n
set.seed(1)

N_approx = length(N_vec)
LIL_n = numeric(N_approx)      # store the true log marginal likelihood
LIL_hat_n = numeric(N_approx)  # store the approximate log marginal likelihood

for (n in 1:N_approx) {
    
    #### generate data
    N = N_vec[n]
    y = rnorm(N, mu, sqrt(sigma_sq))
    
    #### compute posterior parameters
    ybar = mean(y)
    m_n = (N * ybar + w_0 * m_0) / (N + w_0)
    w_n = w_0 + N
    r_n = r_0 + N
    s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2
    
    #### compute true log marginal likelihood for these n data points
    LIL_n[n] = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
        lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
        r_n / 2 * log(s_n)
    
    
    #### begin algorithm
    
    N_iters = 100                 # number of approximations
    def_approx = numeric(N_iters) # storage for default approximations (no r.p.)
    
    
    
    
}



# for n in 1:N { 
#     
#     # (0) generate data, compute: LIL_n
#     # (1) form approximation:     LIL_hat_n
# }









