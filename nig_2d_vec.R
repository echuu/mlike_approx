

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







# automate process for obtaining parameter support

# automate process for calculating c_k

# automate process for calculating closed form integral

# automate product over each dimension for calculating zhat



















