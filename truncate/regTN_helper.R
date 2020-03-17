



# ------------------------------------------------------------------------------
# psi() : - loglikelihood - logprior
# loglikelihood computed using closed form (faster than using dnorm() function)
# logprior computed using TruncatedNormal::dtmvnorm() function
#
psi = function(u, prior) {
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    tau     = prior$tau
    
    beta = u
    
    # psi_u = 0.5 * (N + D) * log(2 * pi * sigmasq)  - 0.5 * D * log(tau) +
    #     1 / (2 * sigmasq) * sum((y - X %*% beta)^2) +
    #     tau / (2 * sigmasq) * sum(beta^2)
    
    loglik = - 0.5 * N * log(2 * pi * sigmasq) - 
        1 / (2 * sigmasq) * sum((y - X %*% beta)^2)
    
    # logTN = TruncatedNormal::dtmvnorm(u, mu = rep(0, D),
    #                  sigma = sigmasq / tau * diag(1, D) , lb = rep(0, D),
    #                  ub = rep(Inf, D), log = T)
    
    logTN = sum(dnorm(u, 0, sqrt(sigmasq / tau), log = T) ) + D * log(2)
    
    psi_u = - loglik - logTN
    
    # matches the calculation below
    # psi_u = -(sum(dnorm(y, X %*% beta, sqrt(sigmasq), log = T)) + 
    #               sum(dnorm(beta, 0, sqrt(sigmasq / tau), log = T)))
    
    
    return(psi_u)
    
} # end of psi() function ------------------------------------------------------




# ------------------------------------------------------------------------------
# lambda() : gradient of psi(u)
# this lambda should be the same as the lambda from previous regression problem
# truncation should not affect any of the computation involved in this
#
lambda = function(u, prior) { 
    
    # note that here we're passing in the posterior params rather than the prior
    # dbeta = Q %*% beta - b
    
    ## TODO: re-write to use the code above, so we can comment out below code
    # verify it's accurate for the normal settings before trying to save
    # on computation
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    tau     = prior$tau
    
    beta = unname(unlist(u))
    # beta = u
    # all quantities below have been calculated before
    # TODO: optimize this part later, can save big-time on computation
    dbeta = 1 / sigmasq * 
        ((t(X) %*% X + tau * diag(D)) %*% beta - t(X) %*% y) 
    
    return(dbeta)
    
} # end of lambda() function ---------------------------------------------------










