

# ------------------------------------------------------------------------------


library(mvtnorm)           # for draws from multivariate normal
library('MCMCpack')        # for rinvgamma() function
library(sn)
library(VGAM)

library(tree)

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)


source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions

source("skew/mv_skew_normal_helper.R")  # load psi(), lambda()


# fixed settings ---------------------------------------------------------------
D = 2
N = 1000 # pseudo-sample size
Omega = diag(1, D)
Sigma = D / N * Omega 
Sigma_inv = solve(Sigma)
alpha = rep(1, D) 
mu_0 = rep(0, D)

prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, 
             alpha = alpha, mu_0 = mu_0)

# -4.3767 for D = 2, N = 500
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) 

# ------------------------------------------------------------------------------

set.seed(1)
J = 10000
N_approx = 10
u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew) # -3.085342




# step through of algorithm to see values of zhat ------------------------------
u_tree = tree(psi_u ~ ., u_df_full)
plot(u_tree)
text(u_tree, cex = 0.8)

plot(u_df_full[,1], u_df_full[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'u1', ylab = 'u2', main = '')
partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


u_df = u_df_full

u_rpart = rpart(psi_u ~ ., u_df)
plot(u_rpart)
text(u_rpart, cex = 0.7)



# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition
u_partition = paramPartition(u_rpart, param_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(u_rpart, u_df, u_partition, D)

n_partitions = nrow(u_partition)     # numebr of partitions 
c_k          = numeric(n_partitions) # constant term for k-th partition
zhat         = numeric(n_partitions) # integral over k-th partition

# (4) compute closed form integral over each partition

lambda_mat = matrix(NA, n_partitions, D)
integral_mat = matrix(NA, n_partitions, D)

for (k in 1:n_partitions) {
    
    # extract "representative point" of the k-th partition
    star_ind = grep("_star", names(param_out))
    u = param_out[k, star_ind] %>% unlist %>% unname
    
    
    # compute lambda_k : gradient of psi, evaluated at u_star
    l_k = lambda(u, prior)       # (D x 1) 
    
    # evaluate e^c_k = e^{psi(u_star)}
    # c_k[k] = exp(-psi(u, prior)) # (1 x 1)
    c_k[k] = exp(-psi(u, prior) + sum(l_k * u)) 
    
    lambda_mat[k,] = l_k
    
    # store each component of the D-dim integral 
    integral_d = numeric(D)      # (D x 1)
    
    for (d in 1:D) {
        
        # find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
            (exp(- l_k[d] * param_out[k, col_id_ub]) - 
                 exp(- l_k[d] * param_out[k, col_id_lb])) 
        
        integral_d[d] = (param_out[k, col_id_ub] - param_out[k, col_id_lb])
        
    } # end of loop computing each of 1-dim integrals
    
    integral_mat[k,] = integral_d
    
    # compute the D-dim integral (product of D 1-dim integrals)
    zhat[k] = prod(c_k[k], integral_d)
    
} # end of for loop over the K partitions


# store the log integral \approx log marginal likelihood
log(sum(zhat))

D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) 

cbind(param_out[,1:4], zhat) %>% cbind(lambda_mat)









# ------------------------------------------------------------------------------
plot(u_samps)
u_df_full = preprocess(u_samps, D)

# skew_rpart = rpart(psi_u ~ ., u_df_full)
# plot(skew_rpart)
# approx_skew

D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) # -4.376731 for D = 2

# ------------------------------------------------------------------------------


# repeat analysis for many D
set.seed(1)
N = 5000 # pseudo-sample size


D_vec = c(2:10)
LIL_D = numeric(length(D_vec))
LIL_D_approx = numeric(length(D_vec))
for (d_i in 1:length(D_vec)) {
    
    D = D_vec[d_i]
    
    Omega = diag(1, D)
    Sigma = D / N * Omega 
    Sigma_inv = solve(Sigma)
    alpha = rep(1, D) 
    mu_0 = rep(0, D)
    
    LIL_D[d_i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5)
    
    J = 5000
    N_approx = 10
    u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame 
    
    u_df_full = preprocess(u_samps, D)
    
    # skew_rpart = rpart(psi_u ~ ., u_df_full)
    # plot(skew_rpart)
    
    approx_skew = approx_lil(N_approx, D, u_df_full, J/N_approx)
    LIL_D_approx[d_i] = mean(approx_skew)
    
}



rbind(LIL_D, LIL_D_approx)


# ------------------------------------------------------------------------------












