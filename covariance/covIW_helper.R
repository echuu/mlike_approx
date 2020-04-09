
## covIW_helper.R

# source("C:/Users/ericc/Dropbox/eric chuu research/GGM/Q1.R")

source("C:/Users/ericc/Dropbox/eric chuu research/GGM/Wishart_InvA_RNG.R")


logmultigamma = function(p, a) {
    f = 0.25 * p * (p - 1) * log(pi)
    for(i in 1:p) { f = f + lgamma( a + 0.5 - 0.5 * i) }
    return(f)
} # end logmultigamma() function -----------------------------------------------



# sampleIW() function
# input: 
#        J     : # of MCMC samples to draw from inverse wishart
#        N     : sample size
#        D_u   : dim of vector u (lower triagular elements of cholesky factor)
#        nu    : prior degrees of freedom
#        S     : sum_n x_n x_n'
#        Omega : prior scale matrix
sampleIW = function(J, N, D_u, nu, S, Omega) {
    
    
    Sigma_post = vector("list", J)  # store posterior samples in matrix form
    L_post     = vector("list", J)  # store lower cholesky factor in matrix form
    post_samps = matrix(0, J, D_u)  # store posterior samples in vector form
    
    for (j in 1:J) {
        
        # Sigma_post_j = solve(Wishart_InvA_RNG(nu + N, S + Omega))
        
        Sigma_post_j = riwish(nu + N, Omega + S)
        
        L_post_j = t(chol(Sigma_post_j)) 
        # store the matrices Sigma, L so that we don't have to reconstruct
        # these matrices later in the functions: 
        # psi(), lambda(), loglik(), cov_loglik()
        
        Sigma_post[[j]] = Sigma_post_j                             # (D x D)
        L_post[[j]]     = L_post_j                                 # (D x D) 
        post_samps[j,]  = L_post_j[lower.tri(L_post_j, diag = T)]  # (D_u x 1)
        
    } # end of sampling loop
    
    
    
    return(list(post_samps = data.frame(post_samps), Sigma_post = Sigma_post, 
                L_post = L_post))
    
    
} # end sampleIW() function ----------------------------------------------------


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

lil = function(param_list) {
    
    N     = param_list$N      # number of observations
    nu    = param_list$nu     # prior degrees of freedom for inv-wishart
    D     = param_list$D      # num cols/rows in covariance matrix Sigma
    Omega = param_list$Omega  # prior scale matrix for inv-wishart
    S     = param_list$S      # sum_n x_n x_n'
    
    logML = logmultigamma(D, (nu + N) / 2) - 0.5 * N * D * log(pi) - 
        logmultigamma(D, nu / 2) + 0.5 * nu * log_det(Omega) - 
        0.5 * (N + nu) * log_det(Omega + S)
    
    return(logML)
    
} # end logML() function -------------------------------------------------------


## maxLogLik() function
maxLogLik = function(Sigma, param_list) {
    
    N     = param_list$N      # number of observations
    D     = param_list$D      # num cols/rows in covariance matrix Sigma
    S     = param_list$S      # sum_n x_n x_n'
    
    loglik = - 0.5 * N * D * log(2 * pi) - 0.5 * N * log_det(Sigma) -
        0.5 * matrix.trace(solve(Sigma) %*% S)
    
    return(loglik)
    
} # end maxLogLik() function ---------------------------------------------------



## cov_loglik() function -------------------------------------------------------
cov_loglik = function(u, params) {
    
    N = params$N
    D = params$D
    S = params$S
    
    L = matrix(0, D, D)             # (D x D) lower triangular matrix
    L[lower.tri(L, diag = T)] = u   # populate lower triangular terms
    
    logDiagL = log(diag(L))         # log of diagonal terms of L
    
    
    # compute loglikelihood using the 'L' instead of 'Sigma' --> allows us
    # to use simplified version of the determinant 
    
    loglik = - 0.5 * N * D * log(2 * pi) - N * sum(logDiagL) - 
        0.5 * matrix.trace(solve(L %*% t(L)) %*% S)
    
    return(loglik)
    
} # end of cov_loglik() function -----------------------------------------------



## cov_logprior() function  ----------------------------------------------------
cov_logprior = function(u, params) {
    
    Omega = params$Omega
    
    nu = params$nu
    D  = params$D  # dimension of covariance matrix Sigma, L
    
    # we do this instead of extracting the covariance matrix explicitly
    # so that we can use the apply function when evaluating psi() for each
    # of the samples 
    
    L = matrix(0, D, D)             # (D x D) lower triangular matrix
    L[lower.tri(L, diag = T)] = u   # populate lower triangular terms
    
    logDiagL = log(diag(L))         # log of diagonal terms of L
    
    # compute log of constant term
    logC = 0.5 * nu * log_det(Omega) - 0.5 * (nu * D) * log(2) - 
        logmultigamma(D, nu / 2)
    
    # compute log of the Jacobian term
    logJacTerm = D * log(2) + sum((D + 1 - 1:D) * logDiagL)
    
    # compute full log prior expression -- note log(det(L)) = sum(logDiagL)
    logprior = logC - (nu + D + 1) * sum(logDiagL) - 
        0.5 * matrix.trace(solve(L %*% t(L)) %*% Omega) +
        logJacTerm
    
    return(logprior)
    
} # end of cov_logprior() function ---------------------------------------------





## psi() function  -------------------------------------------------------------
psi = function(u, params) {
    
    loglik = cov_loglik(u, params)
    logprior = cov_logprior(u, params)
    
    
    return(- loglik - logprior)
    
} # end of psi() function ------------------------------------------------------












# ------------------------------------------------------------------------------
## for now, we just use the constant approximation, can come back and add the 
## gradient if we need to
## TODO: lambda() function




cov_logprior_sigma = function(Sigma, params) {

    # compute the logprior using
    #     (1) sigma  (actual prior on the original model)
    #     (2) L      (induced prior)


    Omega = params$Omega

    nu = params$nu
    D  = params$D  # dimension of covariance matrix Sigma, L


    logC = 0.5 * nu * log_det(Omega) - 0.5 * (nu * D) * log(2) -
        logmultigamma(D, nu / 2)

    logprior = logC - 0.5 * (nu + D + 1) * log_det(Sigma) -
        0.5 * matrix.trace(solve(Sigma) %*% Omega)

    return(logprior)

}
# 
# cov_logprior_L = function(L, params) {
#     
#     # compute the logprior using 
#     #     (1) sigma  (actual prior on the original model)
#     #     (2) L      (induced prior)
#     
#     
#     Omega = params$Omega
#     
#     nu = params$nu
#     D  = params$D  # dimension of covariance matrix Sigma, L
#     
#     
#     # compute log of constant term
#     logC = 0.5 * nu * log_det(Omega) - 0.5 * (nu * D) * log(2) - 
#         logmultigamma(D, nu / 2)
#     
#     # compute log of the Jacobian term
#     logJacTerm = D * log(2) + sum((D + 1 - 1:D) * log(diag(L)))
#     
#     # compute full log prior expression
#     logprior = logC - (nu + D + 1) * log_det(L) - 
#         0.5 * matrix.trace(solve(L %*% t(L)) %*% Omega) +
#         logJacTerm
#     
#     return(logprior)
#     
# }




