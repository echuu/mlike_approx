


library(bridgesampling)



# this function is just -psi(u)
log_density = function(u, data) {
    
    beta = unname(unlist(u[1:p]))
    sigma2 = unname(unlist(u[D]))
    
    sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigmasq), log = T)) + 
        c(a_0 * log(b_0) - p / 2 * log(2 * pi) - 
              0.5 * log_det(V_beta) - lgamma(a_0) -
              (a_0 + p / 2 + 1) * log(sigmasq) - 
              1 / sigmasq * (b_0 + 0.5 * t(beta - mu_beta) %*% V_beta_inv %*% 
                                 (beta - mu_beta)))
}

log_density = function(u, data) {
    -psi(u, data)
}


samples = data.frame(rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), 
                              rep(Inf, D)))
u_df = preprocess(samples, D, prior)

## sample from posterior
u_samp = rtmvnorm(J, c(mu_beta), Q_beta_inv, rep(0, D), rep(Inf, D))
colnames(u_samp) = names(u_df)[1:D]

# prepare bridge_sampler input()
lb = rep(-Inf, D)
ub = rep(Inf, D)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result <- bridgesampling::bridge_sampler(samples = u_samp, log_posterior = log_density,
                                data = prior, lb = lb, ub = ub, silent = TRUE)

# bridge_result <- bridgesampling::bridge_sampler(samples = u_samp, log_posterior = log_density,
#                                                 data = prior, lb = lb, ub = ub, silent = TRUE)


# bridge_result_warp <- bridge_sampler(samples = u_samp, 
#                                      log_posterior = log_density,
#                                      method = 'warp3', data = NULL, 
#                                      lb = lb, ub = ub, silent = TRUE)


# bridge_result$method
# bridge_result$niter
# bridge_result_warp$logml
bridge_result$logml
true_logml


























