


old_logprior = function(theta, params) {
    
    dist = theta[1:d] - beta0
    tau = theta[d+1]
    logPrior = -0.5*(d)*log2Pi + 0.5 * d * log(tau) + 0.5 * log.dettau0 - 
        0.5*tau*t(dist) %*% tau0 %*% dist + 
        dgamma(tau,shape = 0.5*alpha,rate = 0.5*delta,log=TRUE)
    
    return(logPrior)
}



old_loglike = function(theta, params) {
    beta = theta[1:d]
    tau = theta[d+1]
    z = y - X%*%beta
    logLikelihood = -0.5*n*log2Pi + 0.5*n*log(tau) - 0.5*tau*t(z)%*%z
    return(logLikelihood)
}



logpost = function(theta, params) {
    dist = theta[1:d] - beta0
    tau = theta[d + 1]
    logPosterior = -0.5 * (n + d) * log2Pi + 0.5 * (n + d) * log(tau) + 
        0.5 * log.dettau0 -
        0.5 * tau * (t(dist) %*% M %*% dist) - 
        0.5 * tau * c0 + 0.5 * alpha * log(0.5 * delta) - 
        lgamma(0.5 * alpha) + (0.5 * alpha - 1) * log(tau) - 0.5 * delta * tau
    logPosterior
}


old_psi = function(u, params) {
    -old_loglike(u, params) - old_logprior(u, params)
}



old_grad = function(u, l = NULL) {
        beta = u[1:d] %>% unlist %>% unname
        tau = u[d+1]
        
        diff = beta - mu0
        
        g1 = tau * (Lambda0 %*% (beta - mu0) - t(X) %*%(y - X %*% beta))
        g2 = -1/tau * (0.5 * (n + d + alpha) - 1) +
                           0.5 * (delta + sum((y - X%*%beta)^2) +
                                      t(diff) %*% Lambda0 %*% diff)
        return(c(g1, g2))
}



old_hess = function(u, params) {

    beta = u[1:d] %>% unlist %>% unname
    tau = u[d+1]

    h11 = tau * M
    h12 = M %*% beta - XTy - Lambda0 %*% mu0
    h22 = tau^(-2) * (0.5 *(n + d + alpha) - 1)

    H = matrix(0, nrow = d + 1, ncol = d + 1)
    H[1:d, 1:d] = h11
    H[d+1, 1:d] = t(h12)
    H[1:d, d+1] = h12
    H[d+1, d+1] = h22

    return(H)

}


gibbs_radiata = function(Its, BurnIn, fix, initial = NULL,
                         return.log.posterior = FALSE,
                         return.log.likelihood = FALSE) {
    
    # do site by site updates for fair comparison between methods
    T = matrix(nrow = Its - BurnIn, ncol = d + 1)
    
    # inialize from prior
    if (is.null(initial)) {
        tau = rgamma(1, shape = alpha / 2, rate = delta / 2)
        beta = rnorm(d, mean = mu0, sd = sqrt(1 / (tau * diag(tau0))))
    } else{
        beta = initial[1:d]
        tau = intial[d+1]
    }
    
    sh = 0.5 * (n + d + alpha)
    
    sample.vars = which(fix$vars[1:d] == FALSE)
    
    fix.vars = which(fix$vars[1:d] == TRUE)
    if(length(fix.vars) > 0) beta[fix.vars] = fix$values[fix.vars]
    
    if(fix$vars[d + 1] == TRUE) tau = fix$values[d + 1]
    
    sample.tau = !fix$vars[d+1]
    
    for(ItNum in 1:Its){
        
        #visit each parameter in turn
        for(j in sample.vars){
            w = M[j,]%*%(beta-beta0) - M[j,j]*(beta[j]-beta0[j])
            mu = beta0[j] - w/M[j,j]
            sig = sqrt(1/(tau*M[j,j]))
            beta[j] = rnorm(1,mean=mu,sd=sig)
        }
        
        rt = 0.5 * (t(beta-beta0) %*% M %*% (beta - beta0) + c0 + delta)
        if(sample.tau) tau = rgamma(1,shape = sh,rate = rt)
        
        if(ItNum > BurnIn){
            T[ItNum - BurnIn,] = c(beta, tau)
        }
    }		
    
    return(T)
}

