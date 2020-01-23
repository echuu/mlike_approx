


# psi() : negative log posterior
psi = function(u, prior) {
    
    0.5 * t(u) %*% prior$Sigma_inv %*% u
    
}


# lambda() : gradient of psi
lambda = function(u, prior) {
    
    grad(psi, u, prior = prior)
    
}

grad_psi = function(u, prior) {
    
    prior$Sigma_inv %*% u %>% c()
    
}
