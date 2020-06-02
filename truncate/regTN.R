

library(TruncatedNormal)
library(tmg)
library(mvtnorm)

# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# setwd(LEN_PATH)
# 
# source("partition/partition.R")
# source("extractPartition.R")
# source("hybrid_approx.R")

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")     
source("C:/Users/ericc/mlike_approx/truncate/regTN_helper.R")
D = 30
N = 200
I_D = diag(1, D)


# prior mean
mu_0 = rep(0, D)

tau     = 1 / 4          # precision: inverse of variance
sigmasq = 4              # true variance (1 x 1) 

# true value of beta
set.seed(1)
beta = sample(0:10, D, replace = T)
beta = runif(D)
beta = c(runif(D-1, 0, 1), 0)

# generate the regression data -------------------------------------------------

X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
y   = X %*% beta + eps                          # (N x 1) response vector

# compute posterior parameters -------------------------------------------------
Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
Q_beta_inv =  solve(Q_beta)
b          =  1 / sigmasq * t(X) %*% y
mu_beta    =  Q_beta_inv %*% b


# create prior, post objects to be passed into the hml algorithm 
prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
             Q_beta = Q_beta, b = b)
post = list(Q_beta = Q_beta, Q_beta_inv = Q_beta_inv, mu_beta = mu_beta, b = b)


# sample from posterior distribution -------------------------------------------

## TODO: udpate the parameters that are passed into rtmg() function
## i think we need the **precision** matrix, not the covariance matrix

# R       =  rep(0, D)   # linear part R'X
# ff      =  diag(D)     # ?
# gg      =  rep(0, D)   # constraints (truncate to positive orthant)
# initial =  rep(1, D)   # Set initial point for the Markov chain


J        = 3000        # number of MC samples
B        = 10          # number of batch estimates
N_approx = 1           # number of estimates to compute per iteration b

# samples = data.frame(rtmg(J * B,           # number of samples
#                           M = Q_beta_inv,  # (D x D) precision matrix
#                           r = R,           # linear coeff of log-density
#                           initial,         # initial value of markov chain
#                           f = ff,          #
#                           g = gg))         # m-dim with constraints

# plot 2-d samples to verify truncation
# plot(samples)


lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
    0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
    1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)

# lil_0 - log(TruncatedNormal::pmvnorm(rep(0, D), tau / sigmasq * I_D, lb = rep(0, D), ub = rep(Inf, D))[1]) 
(true_logml = lil_0 + D * log(2) + 
    log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                 lb = rep(0, D), ub = rep(Inf, D))[1]))

# run algorithm ----------------------------------------------------------------


set.seed(1)


# TruncatedNormal::dtmvnorm(c(0.8210398, 0.1810452), mu = rep(0, D),
#                  sigma = sigmasq / tau * diag(1, D) , lb = rep(0, D),
#                  ub = rep(Inf, D), log = F)
# 
# prod(dnorm(c(0.8210398, 0.1810452), 0, sqrt(sigmasq / tau), log = F)) * 4


samples = data.frame(rtmvnorm(J * B, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D)))
# samples = data.frame(rtmvnorm(J * B, R, Q_beta_inv, rep(0, D), rep(Inf, D)))
    
u_df = preprocess(samples, D, prior)

# plot(samples)

# u_df %>% head
# source("hybrid_approx_v1.R") # load main algorithm functions
# 
# hml_approx = hml(N_approx, D, u_df, J, prior)

hml_approx = hml_const(1, D, u_df, J, prior)

hml_approx$const_vec

(true_logml = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))


(orig_partition = hml_approx$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
    arrange(desc(perc)) %>% 
    mutate(contrib = logQ_cstar / sum(logQ_cstar)))

K = nrow(orig_partition)

# modified second stage sampling using residuals -------------------------------















# ------------------------------------------------------------------------------

library(profvis)

profvis({
    hml_approx = hml_const(1, D, u_df, J, prior)
})



hml_approx$const_vec # -433.3092 (D = 2), -473.4021 (D = 20)
true_logml # -478.5213

set.seed(1)
profvis({
    reapprox0 = resample_K_simple(hml_approx, K, prior, D)
})

reapprox0 = resample_K_simple(hml_approx, K, prior, D)
log_sum_exp(reapprox0)



# testing the recursive function
n_stages = 2
n_samps = 10
set.seed(1)

start_time <- Sys.time()
reapprox_recursive = resample(hml_approx, prior, n_stages, D, n_samps)
end_time <- Sys.time()
end_time - start_time

log_sum_exp(reapprox_recursive)
length(reapprox_recursive)


# ------------------------------------------------------------------------------

set.seed(1)
start_time <- Sys.time()
reapprox0 = resample_K(hml_approx, K, prior, D)
# log_sum_exp(reapprox0$all_terms)
# length(reapprox0$all_terms)
# K = nrow(orig_partition)
ts_approx_k = vector("list", K) 
ts_approx_terms = vector("list", K) 
for (k in 1:K) {
    
    # print(paste("third stage on partition ", k, sep = ''))
    
    sub_part_k = reapprox0$ss_partitions[[k]]
    
    
    K_sub = nrow(sub_part_k$param_out)
    ts_approx_k[[k]] = resample_K(sub_part_k, K_sub, prior, D, 5)
    
    ts_approx_terms[[k]] = ts_approx_k[[k]]$all_terms
}

end_time <- Sys.time()
end_time - start_time

log_sum_exp(unlist(ts_approx_terms))
length(unlist(ts_approx_terms))

orig_approx = numeric(G)
ss_approx   = numeric(G)
ts_approx   = numeric(G)
for (g in 1:G) {
    
    samples = data.frame(rtmvnorm(J * B, c(mu_beta), Q_beta_inv, 
                                  rep(0, D), rep(Inf, D)))

    u_df = preprocess(samples, D, prior)
    
    hml_approx = hml_const(1, D, u_df, J, prior)
    
    orig_partition = hml_approx$param_out %>%
            dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>%
            dplyr::mutate(perc = n_obs / sum(n_obs)) %>%
            arrange(desc(perc)) %>%
            mutate(contrib = logQ_cstar / sum(logQ_cstar))
    
    K = nrow(orig_partition)
    
    orig_approx[g] = hml_approx$const_vec 
    
    reapprox0 = resample_K(hml_approx, K)
    ss_approx[g] = log_sum_exp(reapprox0$all_terms)
    
    
    ts_approx_k = vector("list", K) 
    ts_approx_terms = vector("list", K) 
    
    for (k in 1:K) {
        
        # print(paste("third stage on partition ", k, sep = ''))
        
        sub_part_k = reapprox0$ss_partitions[[k]]
        
        
        K_sub = nrow(sub_part_k$param_out)
        ts_approx_k[[k]] = resample_K(sub_part_k, K_sub, 5)
        
        ts_approx_terms[[k]] = ts_approx_k[[k]]$all_terms
    }
    
    ts_approx[g] = log_sum_exp(unlist(ts_approx_terms))
    
    
    print(paste('iter ', g, '/', G, ' : ', 
                round(ss_approx[g], 4), ' (err: ', 
                round(abs(true_logml - ss_approx[g]), 4), ', avg: ', 
                round(mean(ss_approx[1:g]), 4), '), ',
                round(ts_approx[g], 4), ' (err: ', 
                round(abs(true_logml - ts_approx[g]), 4), ', avg: ', 
                round(mean(ts_approx[1:g]), 4), '), ', sep = ''))
    
    
        
}

mean(orig_approx) # no re-partioning
mean(ss_approx)   # second stage 
mean(ts_approx)   # third stage
true_logml        # true logML

abs(true_logml - mean(orig_approx))
abs(true_logml - mean(ss_approx))
abs(true_logml - mean(ts_approx))





reapprox$approx
reapprox$subpartitions[[1]]
reapprox$subpartitions[[K]]

hml_approx$const_vec # -473.4021 (D = 20)
reapprox$approx # -480
true_logml # -478.5213

abs(true_logml - hml_approx$const_vec)
abs(true_logml - log_sum_exp(all_terms))

set.seed(1)
G = 50
resample_df = data.frame(min = 1:K, opt = (K-1):0, approx = NA, err = NA)
for (k in 1:K) {
    
    # approx_k = resampleApprox(hml_approx, k)
    
    reapprox = numeric(G)
    for (g in 1:G) {
        reapprox[g] = resampleApprox(hml_approx, k)$approx
        # print(paste('iter ', g, ': ', reapprox[g], sep = ''))
        
    }
    
    approx_k = mean(reapprox)
    
    resample_df[k,]$approx = approx_k
    resample_df[k,]$err = abs(true_logml - approx_k)
    
    print(paste('iter ', k, '/', K, ' : ', round(approx_k, 4), ' (', 
                round(resample_df[k,]$err, 4), ')', sep = ''))
}












