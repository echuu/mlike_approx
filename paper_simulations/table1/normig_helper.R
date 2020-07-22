

psi = function(u, prior) {
    
    y = prior$y
    mu = unname(unlist(u[1]))
    sigmasq = unname(unlist(u[2])) 
    
    r_0 = prior$r_0
    s_0 = prior$s_0
    
    loglik = sum(dnorm(y, mu, sqrt(sigmasq), log = TRUE))
    logprior = dnorm(mu, prior$m_0, sqrt(sigmasq / prior$w_0), log = TRUE) + 
        r_0/2 * log(s_0/2) - lgamma(r_0/2) - (r_0/2 + 1) * log(sigmasq) - 
        s_0/(2*sigmasq)
    # log(dinvgamma(sigma_sq, r_0/2, s_0/2))
        
    return(-loglik-logprior)
}
