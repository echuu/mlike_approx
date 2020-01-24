

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

source("partition/partition.R")  # load partition extraction functions
source("extractPartition.R")
source("hybrid_approx.R")        # load main algorithm functions
source("mvn/mvn_helper.R")       # load psi(), lambda()



D  = 6
N  = 5000
Sigma = diag(1, D)
Sigma = D / N * diag(1, D)
Sigma_inv = solve(Sigma)
# mu_0 = rep(0, D)

prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv)
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500

set.seed(1)
J = 5000
N_approx = 1
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)

approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew)

# plot the tree
library(tree)
u_tree = tree(psi_u ~ ., u_df_full)
plot(u_tree)
text(u_tree, cex = 0.8)

plot(u_df_full[,1], u_df_full[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'u1', ylab = 'u2', main = '')
partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


# integrate psi over (-1, 1) x (-1, 1) using numerical integration
library(pracma)
# n = N
# fun <- function(x, y) exp(-n*x^2*y^4)
D / 2 * log(2 * pi)
fun = function(x, y) exp(-1/2 * (x^2 + y^2))
# fun = function(x) exp(-1/2 * x^2)

result = pracma::integral2(fun, -100, 100, -100, 100, reltol = 1e-50)
log(result$Q) # 1.837877


# need to run through the code in the latter half of mvn.R for code below 
# to run without errors

taylor_int = function(u1, u2) {
    exp(-(lambda1 * (u1 - u_k_star[1]) + lambda2 * (u2 - u_k_star[2])))
}



# (1) use the partitions and numerically compute the integral for each
partition_integral = numeric(n_partitions)
taylor1_integral   = numeric(n_partitions)
taylor1            = numeric(n_partitions)
area_k             = numeric(n_partitions)

taylor2_closed     = numeric(n_partitions)
taylor2_numer      = numeric(n_partitions)

order1_closed      = numeric(n_partitions) # closed form evaluation of integral
order1_numer       = numeric(n_partitions) # numerically integrated

for (k in 1:n_partitions) {
    
    u1_b = param_out[k, 6]  # upper bound of u1
    u1_a = param_out[k, 5]  # lower bound of u1
    
    u2_b = param_out[k, 8]  # upper bound of u2
    u2_a = param_out[k, 7]  # lower bound of u2
    
    # (1) true value of the integral over the k-th partition
    result = integral2(fun, u1_a, u1_b, u2_a, u2_b, reltol = 1e-50)
    partition_integral[k] = result$Q
    
    # (2) compute integral via one term taylor approximation
    u_k_star = param_out[k, star_ind] %>% unlist %>% unname
    
    # e_ck[k] = exp(-psi(u_k_star, prior))
    area_k[k] = (u1_b - u1_a) * (u2_b - u2_a)
    
    taylor1[k] = exp(-psi(u_k_star, prior))
    
    taylor1_integral[k] = taylor1[k] * area_k[k]
    
    
    # (3) compute integral via two term taylor approximation
    
    lambda_k = lambda(u_k_star, prior)
    lambda1  = lambda_k[1]
    lambda2  = lambda_k[2]
    
    
    # (3.1)
    order1_closed[k]  = - 1/lambda1 * (exp(-lambda1 * u1_b) - exp(-lambda1 * u1_a)) * 
        -1/lambda2 * (exp(-lambda2 * u2_b) - exp(-lambda2 * u2_a)) * 
        exp(lambda1 * u_k_star[1]) * exp(lambda2 * u_k_star[2])
    
    order1_closed[k] = - 1/lambda1 * -1/lambda2 * 
        (exp(-lambda1 * u1_b) - exp(-lambda1 * u1_a)) *
        (exp(-lambda2 * u2_b) - exp(-lambda2 * u2_a)) * 
        exp(lambda1 * u_k_star[1]) * exp(lambda2 * u_k_star[2])
    
    taylor2_closed[k] = taylor1[k] * order1_closed[k]
    
    # (3.2)
    order1_numer[k]  = integral2(taylor_int, u1_a, u1_b, u2_a, u2_b, reltol = 1e-50)$Q
    taylor2_numer[k] = taylor1[k] * order1_numer[k]
    
}

log(sum(taylor2_numer))
log(sum(taylor2_closed))


# verify order 1 approximation -- in progress


# check order1_closed, order1_numer
# wolframalpha value for k = 14 --> 0.0000223995

cbind(order1_closed, order1_numer)    # compare closed form integral vs. numer
cbind(taylor2_closed, taylor2_numer)  # compare first order approximations


Z_numer = sum(taylor2_numer)
Z_closed = sum(taylor2_closed)

taylor2_numer
taylor2_closed

log(sum(taylor2_numer))
log(sum(taylor2_closed))

# integral calculated by summing over the K integrals
Z_sum = sum(partition_integral)
log(Z_sum) # 1.836917 --> matches pretty well


log(sum(taylor2_closed))
log(sum(taylor2_numer))



# verify order 0 approximation -- checked!
log(sum(taylor1_integral))
zhat
log(sum(zhat))







# ------------------------------------------------------------------------------


u0 = c(0.5, -0.5) 
u0 = c(runif(1, u1_a, u1_b), runif(1, u2_a, u2_b))
psi(u0, prior)

# taylor approx - 1 term
# psi(u_k, prior)

# taylor approx - 2 term
psi(u_k, prior) + t(grad(psi, u_k, prior = prior)) %*% (u0 - u_k)







