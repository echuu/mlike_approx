

setwd("C:/Users/ericc/mlike_approx/partition")
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

N = 100


#### generate data
y = rnorm(N, mu, sqrt(sigma_sq))

#### compute posterior parameters
ybar = mean(y)
m_n = (N * ybar + w_0 * m_0) / (N + w_0)
w_n = w_0 + N
r_n = r_0 + N
s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2


#### compute true (log) marginal likelihood
LIL = log_p_y = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
    lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
    r_n / 2 * log(s_n)

print(LIL)


log_nig = function(x, sigmasq, mu, lambda, alpha, beta) {
    
    0.5 * (log(lambda) - log(2 * pi * sigmasq)) + alpha * log(beta) - 
        lgamma(alpha) - (alpha + 1) * log(sigmasq) - 
        1 / (2 * sigmasq) * (2 * beta + lambda * (x - mu)^2)
    
}


# ------------------------------------------------------------------------------

## psi_true:     the true negative log posterior, not available in practice, but 
##               we can evaluate it
## psi:          the negative log posterior as described in the notes, 
##               = -loglik - logprior 
## psi_tilde:    approximation of psi as described in the notes
##               = c_k + lambda_k'u


## 1/11 : define generalized version of these functions that take a vector u
#         as input instead of the explicit paramters -- this means that some
#         of these will need to be pre-processed before passed in, but all this
#         will probably require is that they're row-binded (J x D)


#     |  u_11 u_12 ... u_1d  |
#     |   ----- u_2 -----    |
# u = |   ----- u_3 -----    |
#     |   -----  .  -----    |
#     |   -----  .  -----    |
#     |   ----- u_J -----    |


# u = ( mu' , sigma_sq) \in R^D
psi_true_rf = function(u, u_post) {
    
    D = length(u)
    
    mu = u[1:(D-1)]
    sigmasq = u[D]
    
    m_n = u_post$m_n
    w_n = u_post$w_n
    r_n = u_post$r_n
    s_n = u_post$s_n
    
    return(-log_nig(mu, sigmasq, m_n, w_n, r_n / 2, s_n / 2))
}

#### (**) specify psi = -loglik - logprior
psi_rf = function(u, y, prior) {
    
    D = length(u)
    
    mu = u[1:(D-1)]
    sigmasq = u[D]
    
    m_0 = prior$m_0
    w_0 = prior$w_0
    r_0 = prior$r_0
    s_0 = prior$s_0
    
    loglik = sum(dnorm(y, mu, sqrt(sigmasq), log = TRUE))
    logprior = log_nig(mu, sigmasq, m_0, w_0, r_0 / 2, s_0 / 2)
    out = -loglik - logprior
    
    return(out)
}


lambda_rf = function(u, y, prior) {
    
    n = length(y)
    
    mu_star = u[1:(D-1)]
    sigmasq_star = u[D]
    
    
    m_0 = prior$m_0
    w_0 = prior$w_0
    r_0 = prior$r_0
    s_0 = prior$s_0
    
    # in more complicated versions of this problem, these won't be here, since
    # just numerically differentiate, i.e., the code below just turns into
    # grad(psi(u, y, prior)) or something along those lines (above stays same)
    
    lambda1 = -1 / sigmasq_star * sum((y - mu_star)) + 
        w_0 / sigmasq_star * (mu_star - m_0)
    
    lambda2 = n / (2 * sigmasq_star)  - 
        1 / (2 * sigmasq_star^2) * sum((y - mu_star)^2) +
        (r_0 / 2 + 3 / 2) / sigmasq_star - 
        1 / (2 * sigmasq_star^2) * (w_0 * (mu_star - m_0)^2 + s_0)
    
    return(c(lambda1, lambda2))
    
}


# ------------------------------------------------------------------------------


J = 3000  # number of draws from the posterior
D = 2     # dimension of parameter -- in this case u = (mu, sigmasq) \in R^2

# simulate data 

set.seed(1)

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)

u_post = data.frame(mu_post = mu_post, sigma_sq_post = sigma_sq_post)
post = list(m_n = m_n, w_n = w_n, r_n = r_n, s_n = s_n)

psi_u_rf = psi_true_rf(u_post, post)[,1] # get it out of col-vector form

# psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)

# automate process for creating dataframe that partition.R requires for input

# initialize space for u_df : (J x D)
# input for paramPartition() MUST have parameter names u1, u2, ... up
u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u_rf) # (J x 3)

# fit decision tree
### use rpart to fit partition
u_rpart = rpart(psi_u ~ ., u_df)


# obtain the (data-defined) support for each of the parameters
u_support = rbind(c(min(mu_post), max(mu_post)),
                  c(min(sigma_sq_post), max(sigma_sq_post)))

### obtain partition
u_partition = paramPartition(u_rpart, u_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(u_rpart, u_df, u_partition)

n_partitions = nrow(u_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)

prior = list(m_0 = m_0, w_0 = w_0, r_0 = r_0, s_0 = s_0)


for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    # c_k[k] = exp(-psi(param_out[k,]$u1_star, 
    #                   param_out[k,]$u2_star,
    #                   y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    # not ideal to do this here -- fix later (define another function)
    u = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
    
    c_k[k] = exp(-psi_rf(u, y, prior)) # (1 x 1)
    
    l_k = lambda_rf(u, y, prior)
    
    integral_d = numeric(D) # store each component of the D-dim integral 
    
    
    # nothing to refactor in this loop (i think?) since we're just iterating
    # thru each of the integrals and computing an exponential term
    for (d in 1:D) {
        
        # verify these -- these need to be recalculated if the form of param_out
        # changes (if columns get shuffled)
        col_id_lb = 5 + 2 * (d - 1)
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = -1/l_k[d] * 
            exp(-l_k[d] * (param_out[k,col_id_ub] - param_out[k,col_id_lb]))        
        
    }
    
    zhat[k] = prod(c_k[k], integral_d)
    
}

log(sum(zhat))


# testing the wrapper function that iteratively computes the LIL ---------------

approx_lil = function(N_approx, y, y_bar,
                      m_0, w_0, r_0, s_0,
                      m_n, w_n, r_n, s_n) {
    
    #### algorithm: main loop
    
    N_iters = N_approx
    
    # test_out = numeric()
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
        
        u_post = data.frame(mu_post = mu_post, sigma_sq_post = sigma_sq_post)
        post = list(m_n = m_n, w_n = w_n, r_n = r_n, s_n = s_n)
        
        psi_u_rf = psi_true_rf(u_post, post)[,1] # get it out of col-vector form
        
        # psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)
        
        # automate process for creating dataframe that partition.R requires for input
        
        # initialize space for u_df : (J x D)
        # input for paramPartition() MUST have parameter names u1, u2, ... up
        u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u_rf) # (J x 3)
        
        # fit decision tree
        ### use rpart to fit partition
        u_rpart = rpart(psi_u ~ ., u_df)
        
        
        # obtain the (data-defined) support for each of the parameters
        u_support = rbind(c(min(mu_post), max(mu_post)),
                          c(min(sigma_sq_post), max(sigma_sq_post)))
        
        ### obtain partition
        u_partition = paramPartition(u_rpart, u_support)  # partition.R
        
        # organize all data into single data frame --> ready for approximation
        param_out = u_star(u_rpart, u_df, u_partition)
        
        n_partitions = nrow(u_partition)
        c_k = numeric(n_partitions)
        zhat = numeric(n_partitions)
        
        prior = list(m_0 = m_0, w_0 = w_0, r_0 = r_0, s_0 = s_0)
        
        
        for (k in 1:n_partitions) {
            
            # u_star_k = (mu_k, sigma_sq_k)
            # c_k[k] = exp(-psi(param_out[k,]$u1_star, 
            #                   param_out[k,]$u2_star,
            #                   y, m_0, w_0, r_0, s_0)) # (1 x 1)
            
            # not ideal to do this here -- fix later (define another function)
            u = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
            
            c_k[k] = exp(-psi_rf(u, y, prior)) # (1 x 1)
            
            l_k = lambda_rf(u, y, prior)
            
            integral_d = numeric(D) # store each component of the D-dim integral 
            
            
            # nothing to refactor in this loop (i think?) since we're just iterating
            # thru each of the integrals and computing an exponential term
            for (d in 1:D) {
                
                # verify these -- these need to be recalculated if the form of param_out
                # changes (if columns get shuffled)
                col_id_lb = 5 + 2 * (d - 1)
                col_id_ub = col_id_lb + 1
                
                # d-th integral computed in closed form
                integral_d[d] = -1/l_k[d] * 
                    exp(-l_k[d] * (param_out[k,col_id_ub] - param_out[k,col_id_lb]))        
                
            }
            
            zhat[k] = prod(c_k[k], integral_d)
            
        }
        
        
        def_approx[t] = log(sum(zhat))
    }
    
    return(def_approx)
    
} # end of approx_lil()


set.seed(1)
def_approx = approx_lil(100, y, ybar,
                        m_0, w_0, r_0, s_0,
                        m_n, w_n, r_n, s_n) # (100 x 1)


# evaluate the mean, variance of the approximations
mean(def_approx, na.rm = TRUE) # -106.1804
var(def_approx, na.rm = TRUE)  # 4.491305
print(LIL)


# compare approximation with the true value of the log marginal likelihood
approx_default = def_approx
approx_df = data.frame(def = approx_default, lil = LIL) %>% 
    mutate(iter = 1:length(approx_default))

ggplot(approx_df, aes(iter, def)) + geom_point() + 
    geom_hline(aes(yintercept = LIL), 
               linetype = 'dashed', size = 1.5, color = "red")


# ------------------------------------------------------------------------------





# TODO: look into the benefit of doing the following -- as of right now, what we
# have above seems to be sufficient for generalization

# automate process for obtaining parameter support

# automate process for calculating c_k

# automate process for calculating closed form integral

# automate product over each dimension for calculating zhat



















