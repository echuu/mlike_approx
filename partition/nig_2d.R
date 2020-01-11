


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

## note: (1) in harder problems, we only have psi(), psi_tilde()
##       (2) in conjugate examples, we know the posterior in closed form, so we
##           fit the tree using psi_true(), but the approximation part still 
##           uses psi()

#### (*) specify psi_true (in more difficult problems, this one shouldn't be 
#### avaiable, but we're cheating a little to use this when fitting the 
#### decision tree)

psi_true = function(mu, sigma_sq, m_n, w_n, r_n, s_n) {
    
    return(-log_nig(mu, sigma_sq, m_n, w_n, r_n / 2, s_n / 2))
}


#### (**) specify psi = -loglik - logprior
psi = function(mu, sigmasq, y, m_0, w_0, r_0, s_0) {
    
    loglik = sum(dnorm(y, mu, sqrt(sigmasq), log = TRUE))
    logprior = log_nig(mu, sigmasq, m_0, w_0, r_0 / 2, s_0 / 2)
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


# 1/11: testing 1 approximation ------------------------------------------------


J = 3000  # number of draws from the posterior
D = 2     # dimension of parameter -- in this case u = (mu, sigmasq) \in R^2

# simulate data 

set.seed(1)

# (0) sample from mu | sigma_sq, y
mu_post = rnorm(J, m_n, sqrt(sigma_sq / w_n)) # (D x 1)

# (1) sample from sigma_sq | y
sigma_sq_post = MCMCpack::rinvgamma(J, shape = r_n / 2, scale = s_n / 2)


psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)

# automate process for creating dataframe that partition.R requires for input

# initialize space for u_df : (J x D)
# input for paramPartition() MUST have parameter names u1, u2, ... up
u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u) # (J x 3)

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

for (k in 1:n_partitions) {
  
  # u_star_k = (mu_k, sigma_sq_k)
  c_k[k] = exp(-psi(param_out[k,]$u1_star, 
                    param_out[k,]$u2_star,
                    y, m_0, w_0, r_0, s_0)) # (1 x 1)
  
  l_k = lambda(param_out[k,]$u1_star, param_out[k,]$u2_star,
               y, m_0, w_0, r_0, s_0)
  
  integral_d = numeric(D)
  
  for (d in 1:D) {
    
    # verify these
    col_id_lb = 5 + 2 * (d - 1)
    col_id_ub = col_id_lb + 1
    
    integral_d[d] = -1/l_k[d] * 
      exp(-l_k[d] * (param_out[k,col_id_ub] - param_out[k,col_id_lb]))        
    
  }
  
  zhat[k] = prod(c_k[k], integral_d)
  
}

log(sum(zhat))








# ------------------------------------------------------------------------------







#### (***) specify psi_tilde = c_k + lambda_k'u
#### this is calculated within the main loop -- move it out later


## approx_lil() ----------------------------------------------------------------
#### input: 
####         N_approx   : # of marginal likelihood approximations
####         y          : (N x 1) response vector 
####         y_bar      : (1 x 1) mean response
####         m_0        : (D x 1) prior mean
####         w_0        : (1 x 1) scale of covariance
####         r_0        : (1 x 1) shape
####         s_0        : (1 x 1) scale
#### output: def_approx : (N_approx x 1) vector of marginal likelihood approx
####
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
        
        # ***** overflow issues
        psi_u = psi_true(mu_post, sigma_sq_post, m_n, w_n, r_n, s_n)
        
        # input for paramPartition() MUST have parameter names u1, u2, ... up
        u_df = data.frame(u1 = mu_post, u2 = sigma_sq_post, psi_u = psi_u) # (J x 3)
        
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
        for (k in 1:n_partitions) {
            
            # u_star_k = (mu_k, sigma_sq_k)
            #### c_k calculation resulting in very small (1e-48) quantities
            #### that (I think) are causing underflow -> NaN when log(sum(-))
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
            
            # generalize this part for d-th param calulcation
            
            # print(paste("zhat =", zhat[k]))
        }
        
        
        def_approx[t] = log(sum(zhat))
    }
    
    return(def_approx)
    
} # end of approx_lil()
# ------------------------------------------------------------------------------


set.seed(1)
def_approx = approx_lil(100, y, y_bar,
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

# vectorizing u_df initialzation, d'-dim integral calculation













#### log marginal likelihood analysis ------------------------------------------


# (1) for different sample size, compute N_approx LIL approximations
# (2) plot the mean of the approximations over log(sample size)
# (3) fit regression line through these points, verify slope is -d/2

# number of data points for each simulation
N_sample = c(50, 100, 200, 500, 1000, 2000, 10000)

N_sample = c(50, 100, 200, 300, 350)

# number of simulations (number of main loop iterations)
N_sims = length(N_sample)

# number of approximations to compute per iteration of main loop
N_approx = 100


# initialize data storage
LIL_true   = numeric(N_sims)              # true LIL for each sample size
LIL_approx = numeric(N_sims)              # approx LIL for each sample size
LIL_mat    = matrix(0, N_approx, N_sims)  # all approximations stored row-wise

for (t in 1:N_sims) { # main loop
    
    N = N_sample[t]
    
    
    print(paste("sample size =", N, "----------------"))
    
    # generate data
    y = rnorm(N, mu, sqrt(sigma_sq))
    
    #### compute posterior parameters
    ybar = mean(y)
    m_n = (N * ybar + w_0 * m_0) / (N + w_0)
    w_n = w_0 + N
    r_n = r_0 + N
    s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2
    
    
    LIL_true[t] = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
        lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
        r_n / 2 * log(s_n)
    
    LIL_mat[,t] = approx_lil(N_approx, y, y_bar,
                             m_0, w_0, r_0, s_0,
                             m_n, w_n, r_n, s_n) # (100 x 1)
    
    
} # end of main loop


# compute approximate LIL for each sample size by taking column mean of LIL_mat
LIL_approx = colMeans(LIL_mat, na.rm = T)


# plot mean of approximations (for each n) against log(n)
lil_df = data.frame(log_n = N_sample, lil = LIL_approx)

ggplot(lil_df, aes(log_n, lil)) + geom_point()

lm(lil~log_n, lil_df[1:4,])


# plot lil vs log n ------------------------------------------------------------


set.seed(1)

mu = 30
sigma_sq = 4
m_0 = 0
w_0 = 0.05
r_0 = 3
s_0 = 3

N_points = seq(50, 1e3, 100)


#### generate data
LIL_n = numeric(100)
LIL = numeric(length(N_points))

for (n in 1:length(N_points)) {
    
    N = N_points[n]
    
    for (t in 1:100) {
      
      y = rnorm(N, mu, sqrt(sigma_sq))
      
      #### compute posterior parameters
      ybar = mean(y)
      m_n = (N * ybar + w_0 * m_0) / (N + w_0)
      w_n = w_0 + N
      r_n = r_0 + N
      s_n = s_0 + sum((y - ybar)^2) + (N * w_0 / (N + w_0)) * (ybar - m_0)^2
      
      
      #### compute true (log) marginal likelihood
      
      # (LIL = log(p_y)) # -108.877
      
      LIL_n[t] = -(N/2) * log(pi) + 0.5 * (log(w_0) - log(w_n)) + 
        lgamma(r_n / 2) - lgamma(r_0 / 2) + r_0 / 2 * log(s_0) - 
        r_n / 2 * log(s_n)
      
    }
    
    
    LIL[n] = mean(LIL_n)
    
    
}


lil_df = data.frame(log_n = N_points, lil = LIL)

ggplot(lil_df, aes(log_n, lil)) + geom_point()

lm(lil~log_n,lil_df)



# ------------------------------------------------------------------------------



#### ---------------------------------------------------------------------------








