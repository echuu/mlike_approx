


# for now, when moving this into helper function, just load the helper function 
# after loading hybrid_approx.R so that it overrides it
# this function, and psi() both use 'params' arg instead of 'prior' --
# should be switching over to params anyway, it's not always just the prior 
# parameters getting passed into preprocess, psi..
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
    
} # end of preprocess() function




psi = function(u, params) {
    
    # note that the norm(,type = 'F') runs 10x faster than computing the trace
    
    n = params$n # num of rows in X, Y
    p = params$p # num rows of A
    r = params$r # num cols of A, num rows of B^T
    q = params$q # num cols of B^T
    
    sig2 = params$sig2
    del  = params$del
    
    Y = params$Y
    X = params$X
    
    # extract and reshape u to get the matrices: A, B^T
    A_post  = matrix(u[1:(p * r)], p, r)
    Bt_post = matrix(u[(p * r + 1):d], r, q)
    
    const_term = - 0.5 * n * q * log(2 * pi * sig2)
    exp_term   = -1/(2 * sig2) * 
        norm(Y - X %*% A_post %*% Bt_post, type = 'F')^2 + 
        del^2 / (2 * sig2) * 
        (norm(Bt_post, type = 'F')^2 + norm(A_post, type = 'F')^2)
    
    return(const_term + exp_term)
    
} # end of psi() function



### (vectorized) gradient calculation
lambda = function(u, params) {
    
    n = params$n # num of rows in X, Y
    p = params$p # num rows of A
    r = params$r # num cols of A, num rows of B^T
    q = params$q # num cols of B^T
    
    sig2 = params$sig2
    del  = params$del
    
    Y = params$Y
    X = params$X
    
    
    A_post  = matrix(u[1:(p * r)], p, r)
    Bt_post = matrix(u[(p * r + 1):d], r, q)
    B_post  = t(Bt_post)
    
    
    lambda_A  = t(X) %*% Y %*% B_post - 
        t(X) %*% X %*% A_post %*% Bt_post %*% B_post + del^2 * A_post
    
    lambda_Bt = t(A_post) %*% t(X) %*% Y - 
        t(A_post) %*% t(X) %*% X %*% A %*% Bt_post + del^2 * Bt_post
    
    if (nrow(lambda_A)  != nrow(A_post)  | ncol(lambda_A)  != ncol(A_post) |
        nrow(lambda_Bt) != nrow(Bt_post) | ncol(lambda_Bt) != ncol(lambda_Bt)) {
        warning("dimension mismatch in gradient calculation")
    } 
    
    # vectorize the resulting derivatives, column-wise collapse of each matrix
    vec_lambda_A  = c(lambda_A)     # (p*r x 1) vector 
    vec_lambda_Bt = c(lambda_Bt)    # (r*q x 1) vector
    
    return(1 / sig2 * c(vec_lambda_A, vec_lambda_Bt))
    
} # end of lambda() function









