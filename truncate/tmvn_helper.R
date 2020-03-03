


library(numDeriv)

# psi() : negative log posterior
psi = function(u, prior) {
    
    D        =  prior$D
    rho      =  prior$rho
    ld_sigma =  prior$ld_sigma
    
    # 0.5 * t(u) %*% prior$Sigma_inv %*% u
    # prior$D / 2 * log(2 * pi) + 0.5 * log_det(prior$Sigma) +
    #     t(u) %*% prior$Sigma_inv %*% u

    prior$D / 2 * log(2 * pi) + 0.5 * ld_sigma +
        0.5 * t(u) %*% prior$Sigma_inv %*% u
    
    # prior$D / 2 * log(2 * pi) + 0.5 * ld_sigma +
    #     a * sum(u^2) + b * (sum(u^2))^2
} 

# lambda() : gradient of psi -- replace this definition w/ closed
# form definition to speed up the code
lambda = function(u, prior) {
    
    # grad(psi, u, prior = prior)
    prior$Sigma_inv %*% u %>% c()
    
}











