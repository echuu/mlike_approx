



# normal mixture distribution where mixture components are D = 20
# two-dimensional Gaussians, each with 
#     (1) covariance matrix : I_2
#     (2) mean vector       : u_d ~ unif ( [0, 10] x [0, 10] )

D = 2 # dimension of data drawn
n_mix = 20

mu_0 = matrix(c(c(2.18, 5.76), c(8.67, 9.59), c(4.24, 8.48), c(8.41, 1.68), 
                c(3.93, 8.82), c(3.25, 3.47), c(1.70, 0.50), c(4.59, 5.60), 
                c(6.91, 5.81), c(6.87, 5.40), c(5.41, 2.65), c(2.70, 7.88),
                c(4.98, 3.70), c(1.14, 2.39), c(8.33, 9.50), c(4.93, 1.50), 
                c(1.83, 0.09), c(2.26, 0.31), c(5.54, 6.86), c(1.69, 8.11)), 
                nrow = 2)  # mean components of each gaussian 

I_2 = diag(1, D)           # (2 x 2) identity 
w_0 = rep(0.05, 20)        # mixture weights
sigma0 = 0.1               # standard deviation of each gaussian component


prior = list(w = w_0, mu = mu_0, sigma = sigma0, I_D = I_2, n_mix = n_mix)

psi = function(u, prior) {
    
    gauss_vec = numeric(prior$n_mix)
    
    for (j in 1:prior$n_mix) {
        gauss_vec[j] = mvtnorm::dmvnorm(u, mean = prior$mu[,j], 
                                        sigma = prior$sigma * prior$I_D)
    }
        
    return(-log(sum(prior$w * gauss_vec)))
    
}

lambda = function(u, prior) {
    grad(psi, u, prior = prior)
}


# generate data
J = 10000             # number of samples to draw
J_unif = runif(J)
w0_cdf = cumsum(w_0)
mix_id = numeric(J)
for (j in 1:J) {
    mix_id[j] = length(w_0) - sum(J_unif[j] <= w0_cdf) + 1
}

u_samps = matrix(NA, J, 2)
u_samps = data.frame(matrix(NA, J, 2))
names(u_samps) = c("u1", "u2")

for(j in 1:J) {
    u_samps[j,] = mvrnorm(1, mu_0[, mix_id[j]], sigma0[mix_id[j]] * I_2)
}

x11()
ggplot(u_samps, aes(u1, u2)) + geom_point(size = 0.8)

u_df = preprocess(u_samps, D, prior) # takes long time to run -- maybe write in STAN

library(tree)
gm_tree = tree(psi_u ~ ., u_df)
partition.tree(gm_tree)

# run algorithm
N_approx = 10
J = 1e4
approx_out = hybrid_mlik(N_approx, D, u_df, J / N_approx, prior) 


# plot partition
x11()
plotPartition(u_df, approx_out$verbose_partition)

# examine approximations
approx_out$hybrid_vec %>% mean
approx_out$taylor_vec %>% mean
approx_out$const_vec  %>% mean


# ------------------------------------------------------------------------------

library(pracma)

fun = function(x1, x2) {
    
    # gauss_vec = numeric(prior$n_mix)
    # 
    # for (j in 1:prior$n_mix) {
    #     gauss_vec[j] = mvtnorm::dmvnorm(c(x1, x2), mean = c(prior$mu[,j]), 
    #                                     sigma = prior$sigma * prior$I_D)
    #     
    #     # gauss_vec[j] = dnorm(x1, mean = prior$mu[1,j], sd = prior$sigma) * 
    #     #     dnorm(x2, mean = prior$mu[2,j], sd = prior$sigma)
    #     
    # }
    
    return(sum(prior$w * 
                   apply(prior$mu, 2, dmvnorm, x = c(x1, x2), 
                         sigma = prior$sigma * prior$I_D)))
}
result = integral2(fun, 0, 10, -1, 11, reltol = 1e-50)
result

log(result$Q) # -1.223014 for n = 1000





# ------------------------------------------------------------------------------


# true value of the normalizing constant ??
install.packages("bridgesampling")
library(bridgesampling)


jags_H1 <- jags(data = list(u_samps),
                parameters.to.save = c("u1", "u2"),
                model.file = textConnection(code_H1), n.chains = 3,
                n.iter = 16000, n.burnin = 1000, n.thin = 1)

log_posterior_H1 = function(pars, data) {
    
    gauss_vec = numeric(pars$n_mix)
    
    for (j in 1:pars$n_mix) {
        gauss_vec[j] = mvtnorm::dmvnorm(data$u, mean = pars$mu[,j], 
                                        sigma = pars$sigma * pars$I_D)
    }
    
    return(log(sum(prior$w * gauss_vec)))
    
}
lb_H1 = c(-Inf, Inf)
ub_H1 = c(-Inf, Inf)
bridge_H1 = bridge_sampler(samples = jags_H1,
                           log_posterior = log_posterior_H1,
                           data = list(u_samps),
                           lb = lb_H1, ub = ub_H1)

print(bridge_H1)



