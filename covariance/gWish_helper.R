

sampleGW = function(J, D_u, D_0, G, b, N, V, S, edgeIndex) {
    
    Omega_post    = vector("list", J) # store posterior samples in matrix form
    Lt_post       = vector("list", J) # store lower cholesky factor
    post_samps_0  = matrix(0, J, D_0) # store ENTIRE upper diag in vector form
    post_samps    = matrix(0, J, D_u) # store NONZERO upper diag in vector form
    
    
    
    Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
    
    Lt_post_j = chol(Omega_post)
    
    for (j in 1:J) {
        
        ## we store the next 2 quantities purely for testing so we can ---------
        ## verify that we can reconstruct these two matrices
        
        Lt_post_j = chol(Omega_post[,,j])
        Lt_post[[j]] = Lt_post_j
        
        # collapse upper triangular matrix into a vector
        Lt_upper_vec = Lt_post_j[upper.tri(Lt_post_j, diag = T)]    # (D_0 x 1)
        # store entire vector (includes 0 elements)
        post_samps_0[j,]  = Lt_upper_vec                            # (D_0 x 1)
        
        ##
        ## ---------------------------------------------------------------------

        # this is what we observe in the first D_u columns of u_df later
        # store nonzero vector (this is used to fit the tree)       # (D_u x 1)
        post_samps[j,] = Lt_upper_vec[edgeIndex]
    }
    
    
    
    
}


























