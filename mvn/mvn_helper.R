


log_det = function(xmat) {
    return(c(determinant(xmat, logarithm = T)$modulus))
}



# psi() : negative log posterior
psi = function(u, prior) {
    
    0.5 * t(u) %*% prior$Sigma_inv %*% u %>% c
    
}


# lambda() : gradient of psi
lambda = function(u, prior) {
    
    
    grad(psi, u, prior = prior)
    
    
}

