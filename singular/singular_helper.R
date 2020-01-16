


# psi() : negative log posterior
psi = function(u, N) {
    
    return(N * u[1]^2 * u[2]^4)
    
}


# lambda() : gradient of psi
lambda = function(u, N) {
    
    l_1 = 2 * N * u[1] * u[2]^4
    l_2 = 4 * N * u[1]^2 * u[2]^3
    
    return(c(l_1, l_2))
}











