

# psi() : negative log posterior
psi = function(u, prior) {
    N = prior$N
    return(N * u[1]^2 * u[2]^4 * u[3]^2 * u[4]^2)
}



