

process_samps = function(Omega){
    Lt = chol(Omega)
    Lt_vec = Lt[upper.tri(Lt, diag = T)]
    return(Lt_vec)
}


sampleGW = function(J, D_0, G, b, N, V, S) {
    
    Omega_post    = vector("list", J) # store posterior samples in matrix form
    # Lt_post       = vector("list", J) # store lower cholesky factor
    # post_samps_0  = matrix(0, J, D_0) # store ENTIRE upper diag in vector form
    # post_samps    = matrix(0, J, D_u) # store NONZERO upper diag in vector form
    
    Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
    
    # Lt_post_j = chol(Omega_post)
    
    # for (j in 1:J) {
    #     
    #     ## we store the next 2 quantities purely for testing so we can ---------
    #     ## verify that we can reconstruct these two matrices
    #     
    #     Lt_post_j = chol(Omega_post[,,j])
    #     Lt_post[[j]] = Lt_post_j
    #     
    #     # collapse upper triangular matrix into a vector
    #     Lt_upper_vec = Lt_post_j[upper.tri(Lt_post_j, diag = T)]    # (D_0 x 1)
    #     # store entire vector (includes 0 elements)
    #     post_samps_0[j,]  = Lt_upper_vec                            # (D_0 x 1)
    # }
    
    post_samps = t(apply(Omega_post, 3, process_samps))
    # return(list(samp_array = Omega_post, post_samps = post_samps_0))
    return(post_samps)
}




gwish_loglik = function(u, params) {
    
    N   = params$N
    D   = params$D
    S   = params$S
    b   = params$b       # degree of freedom for G-wishart distribution
    V   = params$V # scale matrix for G-wishart distribution
    
    Lt = matrix(0, D, D)              # (D x D) lower triangular matrix
    Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
    
    # recall: Sigma^(-1) = LL' 
    loglik = - 0.5 * N * D * log(2 * pi) + N * log_det(Lt) - 
        0.5 * matrix.trace(t(Lt) %*% Lt %*% S)
    
    return(loglik)
}

gwish_logprior = function(u, params) {
    
    G   = params$G
    D   = params$D
    b   = params$b       # degree of freedom for G-wishart distribution
    V   = params$V # scale matrix for G-wishart distribution
    
    Lt = matrix(0, D, D)              # (D x D) lower triangular matrix
    Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
    
    logprior = - gnorm(G, b, V, 100) + (b - 2) * sum(log(diag(Lt))) - 
        0.5 * matrix.trace(t(Lt) %*% Lt %*% V)
    
    return(logprior)
}



## psi() function  -------------------------------------------------------------
psi = function(u, params) {
    
    loglik = gwish_loglik(u, params)
    logprior = gwish_logprior(u, params)
    
    
    return(- loglik - logprior)
    
} # end of psi() function ------------------------------------------------------


























