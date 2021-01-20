

# source("C:/Users/ericc/Dropbox/eric chuu research/GGM/HIWsim.R")


## sampleHIW() function
#  J       : number of posterior samples to generate
#  D_u     : dim of the sample u that is fed into the tree (# nonzero in L')
#  D_0     : num entries on diagonal and upper diagonal (can remove later)
#  b       : df for HIW prior
#  N       : sample size
#  V       : (D x D) p.d. matrix for HIW Prior
#  S       : Y'Y, a (D x D) matrix
sampleHIW = function(J, D_u, D_0, G, b, N, V, S, edgeIndex) {
    
    Sigma_post    = vector("list", J) # store posterior samples in matrix form
    Omega_post    = vector("list", J) # store posterior samples in matrix form
    Lt_post       = vector("list", J) # store lower cholesky factor
    
    post_samps_0  = matrix(0, J, D_0) # store ENTIRE upper diag in vector form
    post_samps    = matrix(0, J, D_u) # store NONZERO upper diag in vector form
    
    for (j in 1:J) {
        
        # perform one draw from HIW(b + N, V + S)
        HIW_draw = HIWsim(G, b + N, V + S)
        draw_Sigma = HIW_draw$Sigma
        draw_Omega = HIW_draw$Omega
        
        Lt_post_j = chol(draw_Omega) # upper cholesky factor for Omega = LL'
        
        Sigma_post[[j]] = draw_Sigma                                # (D x D)
        Omega_post[[j]] = draw_Omega                                # (D x D)
        Lt_post[[j]]    = Lt_post_j                                 # (D x D) 
        
        # collapse upper triangular matrix into a vector
        Lt_upper_vec = Lt_post_j[upper.tri(Lt_post_j, diag = T)]    # (D_0 x 1)
        
        # store entire vector (includes 0 elements)
        post_samps_0[j,]  = Lt_upper_vec                            # (D_0 x 1)
        
        # store nonzero vector (this is used to fit the tree)       # (D_u x 1)
        post_samps[j,] = Lt_upper_vec[edgeIndex]
        
    } # end of sampling loop
    
    
    ## note: Both post_samps_0 and post_samps can be used to re-construct
    ##       Lt_post. The former just needs to be stuck into the upper.tri 
    ##       elements of a (D x D) matrix, while the latter requires the 
    ##       edgeIndex logical vector to index into the correct elements
    
    return(list(post_samps   = data.frame(post_samps), 
                post_samps_0 = data.frame(post_samps_0),
                Lt_post      = Lt_post,
                Sigma_post   = Sigma_post, 
                Omega_post   = Omega_post, Lt_post = Lt_post))
    
} # end sampleHIW() function ---------------------------------------------------



## preprocess() function
#
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





## maxLogLik() function --------------------------------------------------------
# maximized log-likelihood, i.e., log-likelihood evaluated at the true Omega
maxLogLik = function(Omega, params) {

    N     = params$N      # number of observations
    D     = params$D      # num cols/rows in covariance matrix Sigma
    S     = params$S      # sum_n x_n x_n'
    
    loglik = - 0.5 * N * D * log(2 * pi) + 0.5 * N * log_det(Omega) -
        0.5 * matrix.trace(Omega %*% S)
    
    return(loglik)
        
} # end maxLogLik() function ---------------------------------------------------





# ## HIW_loglik() function -----------------------------------------------------
# HIW_loglik_old = function(u, params) {
#     
#     N = params$N
#     D = params$D
#     S = params$S
#     
#     Lt = matrix(0, D, D)              # (D x D) lower triangular matrix
#     Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
#     
#     # recall: Sigma^(-1) = LL'
# 
#     logprior = - 0.5 * N * D * log(2 * pi) + N * log_det(Lt) - 
#         0.5 * matrix.trace(t(Lt) %*% Lt %*% S)
#     
#     return(logprior)
#     
# } # end HIW_loglik() function ------------------------------------------------





HIW_loglik = function(u, params) {
    
    N   = params$N
    D   = params$D
    D_0 = params$D_0
    S   = params$S
    
    Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
    Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
    Lt_vec_0[edgeInd] = u
    Lt[upper.tri(Lt, diag = T)] = Lt_vec_0   # populate lower triangular terms
    
    # recall: Sigma^(-1) = LL'
    
    logprior = - 0.5 * N * D * log(2 * pi) + N * log_det(Lt) - 
        0.5 * matrix.trace(t(Lt) %*% Lt %*% S)
    
    return(logprior)
    
} # end HIW_loglik() function --------------------------------------------------


## HIW_logprior() function -----------------------------------------------------
HIW_logprior = function(u, params) {
    
    # steps:
    # (1) extract upper diagonal entries
    # (2) extract diagonal entries
    # (3) compute nu = (nu_1,...,nu_p)
    # (4) compute upper diagonal part of log prior
    # (5) compute diagonal part of log prior
    
    D   = params$D              # number of rows/cols in Sigma/Omega
    D_0 = params$D_0            # num entries on diagonal and upper diagonal
    b   = params$b              # degrees of freedom
    
    edgeInd  = params$edgeInd   # indicator for present edges in the graph F
    upperInd = params$upperInd  # indicator for upper diagonal edges
    
    Lt = matrix(0, D, D)
    Lt_vec_0 = numeric(D_0)
    Lt_vec_0[edgeInd] = u
    Lt[upper.tri(Lt, diag = T)] = Lt_vec_0 # reconstruct Lt (upper tri matrix)
    
    #### (1) upper diagonal entries
    u_upper_all = Lt[upper.tri(Lt)] # extract upper diagonal entries
    u_upper = u_upper_all[upperInd] # keep only entries that have edge in G
    
    #### (2) diagonal entries
    u_diag = diag(Lt)               # extract diagonal entries, all included
    
    #### (3) compute nu_i, i = 1,..., d
    # compute nu_i (i = 1,...,D) by counting nonzero elements 
    # in each row of Lt - 1
    # recall: the i-th row of Lt has exactly nu_i + 1 nonzero entries
    nu = rowSums(Lt != 0) - 1 
    nu[D] = 0
    
    #### (4) compute upper diagonal part of log prior
    upper_diag_prior = sum(- 0.5 * log(2 * pi) - 0.5 * u_upper^2)
    
    #### (5) compute diagonal part of log prior
    diag_prior = sum(- 0.5 * (b + nu) * log(2) - lgamma(0.5 * (b + nu)) + 
                         (b + nu - 2) * log(u_diag) - 
                         0.5 * u_diag^2 + log(2 * u_diag))
    
    # log prior, as shown in (8) of 3.2.2 HW Induced Cholesky Factor Density
    logprior = upper_diag_prior + diag_prior
    
    return(logprior)
    
    
} # end HIW_logprior() function ------------------------------------------------





## psi() function  -------------------------------------------------------------
psi = function(u, params) {
    
    loglik = HIW_loglik(u, params)
    logprior = HIW_logprior(u, params)
    
    
    return(- loglik - logprior)
    
} # end of psi() function ------------------------------------------------------





# end HIW_helper.R 
