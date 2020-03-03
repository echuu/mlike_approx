


## TODO: lil() : compute the TRUE log marginal likelihood
lil = function(prior, post) {
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    lambda  = prior$lambda
    
    b          = post$b
    mu_beta    = post$mu_beta
    Q_beta     = post$Q_beta
    
    log_py = 0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
        0.5 * D * log(lambda) - 0.5 * log_det(Q_beta) - 
        1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
    
    return(log_py)
    
} # end of lil() function


## TODO: psi() : compute (- loglike - logprior)
## verify this expression with R functions using dnorm(, log = T)
psi = function(u, prior) {
    
    y = prior$y
    X = prior$X
    
    N = prior$N
    D = prior$D
    
    sigmasq = prior$sigmasq
    lambda  = prior$lambda
    
    beta = u
    
    psi_u = 0.5 * (N + D) * log(2 * pi * sigmasq)  - 0.5 * D * log(lambda) + 
        1 / (2 * sigmasq) * sum((y - X %*% beta)^2) + 
        lambda / (2 * sigmasq) * sum(beta^2)
    
    return(psi_u)
    
} # end of psi() function


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
    lambda  = prior$lambda
    
    # beta = unname(unlist(u))
    beta = u
    # all quantities below have been calculated before
    # TODO: optimize this part later, can save big-time on computation
    dbeta = 1 / sigmasq * 
        ((t(X) %*% X + lambda * diag(D)) %*% beta - t(X) %*% y) 
    
    return(dbeta)
        
} # end of lambda() functino




## mvn_ig_helper.R file has a preprocess() that overwrites the generic preprocess()
## function -- don't think that needs to be done here, since just going to be 
## using R to sample from the normal (posterior) distribution
















