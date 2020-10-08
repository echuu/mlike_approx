


loglik = function(beta) {
    - 0.5 * N * log(2 * pi * sigmasq) - 
        1 / (2 * sigmasq) * sum((y - X %*% beta)^2)
}


logZ = -9999 # initial value for log marignal likelihood

J = 100

# sample J points from the prior
pr_samp = rtmvnorm(J, rep(0, D), 1 / tau * sigmasq * diag(1,D), 
                   rep(0, D), rep(Inf, D))

j = 0
H = 0     # information, initially 0
end = 2   # termination condition: nest = end * J * H

# begin nested sampling loop ---------------------------------------------------
set.seed(1)
while (j <= end * J * H) {
    print(paste("iter ", j, ": logZ = ", logZ, sep = ''))
    # (0) outermost interval of prior mass
    logw = log(1 - exp(-1 / J))
    
    # (1) find worst sample (lowest likelihood)
    loglik_vec = apply(pr_samp, 1, loglik) # evaluate loglik for all prior draws
    worst = which.min(loglik_vec)          # index of worst prior
    log_wt = logw + loglik_vec[worst]
    
    # (2) update evidence Z and information H
    logZ_new = log_sum_exp(c(logZ, log_wt))
    H = exp(log_wt - logZ_new) * loglik_vec[worst] + 
        exp(logZ - logZ_new) * (H + logZ) - logZ_new
    logZ = logZ_new
    
    # kill worst object in favor of copy of different survive
    # replace point of lowest likelihood by new one drawn from within
    # L(theta) > L_i, in proportion to the prior pi(theta)
    
    logL_star = loglik_vec[worst] # new likelihood constraint
    logL_new = logL_star
    # draw from prior until we find a prior to over-write the worst prior
    while (logL_new <= logL_star) {
        # print(1)
        theta = rtmvnorm(1, rep(0, D), 1 / tau * sigmasq * diag(1,D), 
                         rep(0, D), rep(Inf, D))
        logL_new = loglik(theta)
    }
    
    # replace worst theta with the new theta
    pr_samp[worst,] = theta
    loglik_vec[worst] = logL_new
    
    # shrink interval
    logw = logw - 1 / J
    
    j = j + 1
} # end of nested sampling loop ------------------------------------------------

logw = -j / J - log(J)
for (i in 1:N) {
    
    logZnew = log_sum_exp(logZ, )
    
    
    
}

LIL


