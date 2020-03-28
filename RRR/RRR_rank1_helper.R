





loglike_true = function(a_0, b_0, params) {
    
    p = params$p  # dimension of a
    q = params$q  # dimension of b
    
    sig2 = params$sig2
    Y = params$Y
    
    - 0.5 * p * q * log(2 * pi * sig2) - 
        1 / (2 * sig2) * norm(Y - a_0 %*% t(b_0), type = 'F')^2
    
} # end of loglike_true() function ---------------------------------------------



rrr_loglik = function(u, params) {
    
    p = params$p  # dimension of a
    q = params$q  # dimension of b
    
    sig2 = params$sig2
    Y = params$Y
    
    a = u[1:p]      # extract a (p x 1) from the posterior sample u
    b = tail(u, q)  # extract b (q x 1) from the posterior sample u
    
    - 0.5 * p * q * log(2 * pi * sig2) - 
        1 / (2 * sig2) * norm(Y - a %*% t(b), type = 'F')^2
    
    
} # end of rrr_loglik() function -----------------------------------------------





rrr_logprior = function(u, params) {
    
    p = params$p  # dimension of a
    q = params$q  # dimension of b
    
    sig2 = params$sig2
    del  = params$del
    
    a = u[1:p]      # extract a (p x 1) from the posterior sample u
    b = tail(u, q)  # extract b (q x 1) from the posterior sample u
    
    - 0.5 * (p + q) * log(2 * pi * sig2) + 0.5 * (p + q) * log(del^2) - 
        del^2 / (2 * sig2) * (sum(a^2) + sum(b^2))
    
} # end of rrr_logprior() function ---------------------------------------------




psi = function(u, params) {
    
    p = params$p  # dimension of a
    q = params$q  # dimension of b
    
    loglik = rrr_loglik(u, params)
    logprior = rrr_logprior(u, params)
    
    return(-loglik - logprior)
    
} # end of psi() function ------------------------------------------------------




lambda = function(u, params) {
    
    p = params$p  # dimension of a
    q = params$q  # dimension of b

    sig2 = params$sig2
    del  = params$del
    
    Y = params$Y
    
    a = u[1:p]      # extract a (p x 1) from the posterior sample u
    b = tail(u, q)  # extract b (q x 1) from the posterior sample u
    
    ata = sum(a^2)  # a'a (1 x 1) scalar
    btb = sum(b^2)  # b'b (1 x 1) scalar
    
    

    lambda_a = btb * a - Y %*% b + del^2 * a  # (p x 1)
    lambda_b = ata * b - Y %*% a + del^2 * b  # (q x 1)
    
    
    # return vectorized version of the gradient
    return(1 / sig2 * c(lambda_a, lambda_b))
    
} # end of lambda() function











