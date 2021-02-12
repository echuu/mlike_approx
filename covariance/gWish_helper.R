

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
    post_samps = t(apply(Omega_post, 3, process_samps))
    return(post_samps)
}



preprocess = function(post_samps, D, params) {
    
    psi_u = apply(post_samps, 1, psi, params = params) %>% unname() # (J x 1)
    
    # (1.2) name columns so that values can be extracted by partition.R
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(post_samps, psi_u) # J x (D + 1)
    names(u_df) = u_df_names
    
    
    return(u_df)
    
} # end of preprocess() function -----------------------------------------------




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
    b   = params$b # degree of freedom for G-wishart distribution
    V   = params$V # scale matrix for G-wishart distribution
    
    Lt = matrix(0, D, D)              # (D x D) lower triangular matrix
    Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
    
    logprior = - gnorm(G, b, V, 100) + (b - 2) * sum(log(diag(Lt))) - 
        0.5 * matrix.trace(t(Lt) %*% Lt %*% V)
    
    return(logprior)
}


gwish_logprior = function(u, params) {
    
    D   = params$D
    b   = params$b # degree of freedom for G-wishart distribution
    V   = params$V # scale matrix for G-wishart distribution
    nu  = params$nu
    
    Lt = matrix(0, D, D)              # (D x D) upper triangular matrix
    Lt[upper.tri(Lt, diag = T)] = u   # populate upper triangular terms
    
    x = (b - 2) * sum(log(diag(Lt))) - 0.5 * matrix.trace(t(Lt) %*% Lt %*% V) + 
        D * log(2) + sum((nu + 1) * log(diag(Lt)))
    
    sum((b + nu - 1) * log(diag(Lt))) - 0.5 * matrix.trace(t(Lt) %*% Lt %*% V) + 
        D * log(2)
    
}


## psi() function  -------------------------------------------------------------
slow_psi = function(u, params) {
    
    # loglik = gwish_loglik(u, params)
    logprior = gwish_logprior(u, params)
    
    return(- logprior)
    
} # end of psi() function ------------------------------------------------------


























