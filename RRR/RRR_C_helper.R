





### preprocess() ---------------------------------------------------------------
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



### loglik_true() --------------------------------------------------------------
## evaluate the loglikelihood of the following regression model:
##     Y ~ MN (XAB^T, I_n, sig2 * I_q)
##     vec(y) ~ N(vec(XAB^T), sig2 kron(I_q, I_n))
## 
loglik_true = function(A_0, B_0, params) {
    
    n = params$n
    q = params$q
    
    Y = params$Y
    X = params$X
    
    sig2 = params$sig2
    
    loglik = (-n * q / 2) * log(2 * pi * sig2) -
        1 / (2 * sig2) * norm(Y - X %*% A_0 %*% t(B_0), type = 'F')^2
    
    return(loglik)
    
} # end loglik_true function ---------------------------------------------------



### loglik() -------------------------------------------------------------------
## evaluate the loglikelihood of the following regression model:
##     Y ~ MN (XAB^T, I_n, sig2 * I_q)
##     vec(y) ~ N(vec(XAB^T), sig2 kron(I_q, I_n))
## 
rrr_loglik = function(u, params) {
    
    p = params$p
    q = params$q
    
    n = params$n
    r = params$r
    
    d = params$d
    
    Y = params$Y
    X = params$X
    
    sig2 = params$sig2
    
    # likelihood evaluated using the posterior samples
    A_post  = matrix(u[1:(p * r)], p, r)
    Bt_post = matrix(u[(p * r + 1):d], r, q)
    
    C_post = matrix(u, p, q)
    
    loglik = (-n * q / 2) * log(2 * pi * sig2) -
        1 / (2 * sig2) * norm(Y - X %*% C_post, type = 'F')^2
    
    return(loglik)
} # end rrr_loglik function ----------------------------------------------------



### lambda() -- (vectorized) gradient calculation ------------------------------
## input : 
## u : (d x 1) vector, d = r * p + r * q
## output: 
## lambda(u) : (d x 1) vector -- vectorized version of the actual gradient
##
psi = function(u, params) {
    
    # note that the norm(,type = 'F') runs 10x faster than computing the trace
    
    n = params$n # num of rows in X, Y
    p = params$p # num rows of A
    r = params$r # num cols of A, num rows of B^T
    q = params$q # num cols of B^T
    d = params$d
    
    sig2 = params$sig2
    del  = params$del
    
    Y = params$Y
    X = params$X
    
    loglik = rrr_loglik(u, params)
    # logprior = rrr_logprior(u, params)
    
    return(-loglik)
    
    # return(const_term + exp_term)
    
} # end of psi() function



### lambda() -- (vectorized) gradient calculation ------------------------------
## input : 
## u : (d x 1) vector, d = r * p + r * q
## output: 
## lambda(u) : (d x 1) vector -- vectorized version of the actual gradient
##
lambda = function(u, params) {
    
    p = params$p # num rows of A
    q = params$q # num cols of B^T
    
    sig2 = params$sig2
    
    XtX = params$XtX
    Xty = params$Xty
    
    C = matrix(u, p, q)
    
    lambda_C = 1 / sig2 * (XtX %*% C - Xty)
    
    
    # return vectorized version of the gradient
    return(c(lambda_C))
    
} # end of lambda() function











