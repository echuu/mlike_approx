


## lil() : compute the TRUE log marginal likelihood
lil = function(prior, post) {
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    tau     = prior$tau
    
    b          = post$b
    mu_beta    = post$mu_beta
    Q_beta     = post$Q_beta
    
    log_py = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
        0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
        1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
    
    return(log_py)
    
} # end of lil() function


## psi() : compute (- loglike - logprior)
## 
psi = function(u, prior) {
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    tau     = prior$tau
    
    beta = u
    
    psi_u = 0.5 * (N + D) * log(2 * pi * sigmasq)  - 0.5 * D * log(tau) +
        1 / (2 * sigmasq) * sum((y - X %*% beta)^2) +
        tau / (2 * sigmasq) * sum(beta^2)
    
    # matches the calculation below
    # psi_u = -(sum(dnorm(y, X %*% beta, sqrt(sigmasq), log = T)) + 
    #               sum(dnorm(beta, 0, sqrt(sigmasq / tau), log = T)))
    
    
    return(psi_u)
    
} # end of psi() function


psi1 = function(u, prior) {
 
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    tau     = prior$tau
    
    beta = u
    
    
    
     
    return(0)
}


## TODO: lambda() : compute gradient of psi()
## verify this expression matches with grad() function output
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
        
} # end of lambda() function




## mvn_ig_helper.R file has a preprocess() that overwrites the generic preprocess()
## function -- don't think that needs to be done here, since just going to be 
## using R to sample from the normal (posterior) distribution
















