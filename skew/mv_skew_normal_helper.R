

psi = function(u, prior) {
    
    0.5 * t(u) %*% prior$Sigma_inv %*% u - 
        log(0.5 * (1 + erf(sum(prior$alpha * u) / sqrt(2))))
    
    
}

lambda =  function(u, prior) {
    grad(psi, u, prior = prior)
}

# ------------------------------------------------------------------------------




# ----------------------------------------------------------------------


#psi_skew = function(u) {

#    0.5 * t(u) %*% Sigma_inv %*% u - 
#        log(0.5 * (1 + erf(sum(alpha * u) / sqrt(2))))

#}

#lambda =  function(u) {
#    grad(psi_skew, u)
#}








