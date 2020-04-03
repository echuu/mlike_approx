
## covIW_helper.R

source("C:/Users/ericc/Dropbox/eric chuu research/GGM/Q1.R")


logmultigamma = function(p, a) {
    f = 0.25 * p * (p - 1) * log(pi)
    for(i in 1:p) { f = f + lgamma( a + 0.5 - 0.5 * i) }
    return(f)
} # end logmultigamma() function -----------------------------------------------


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


## TODO: maxLogLik() function
maxLogLik = function(Sigma, param_list) {
    
    N     = param_list$N      # number of observations
    D     = param_list$D      # num cols/rows in covariance matrix Sigma
    S     = param_list$S      # sum_n x_n x_n'
    
    loglik = - 0.5 * N * D * log(2 * pi) - 0.5 * N * log_det(Sigma) -
        0.5 * matrix.trace(solve(Sigma) %*% S)
    
    return(loglik)
    
} # end maxLogLik() function ---------------------------------------------------






## TODO: log_prior() function



## TODO: cov_loglik() function



## TODO: psi() function


















# ------------------------------------------------------------------------------
## for now, we just use the constant approximation, can come back and add the 
## gradient if we need to
## TODO: lambda() function




