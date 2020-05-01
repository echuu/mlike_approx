
# setwd("C:/Users/ericc/Dropbox/eric chuu research/GGM/")

library(BDgraph)
source("C:/Users/ericc/Dropbox/eric chuu research/GGM/makedecompgraph.R")
source("C:/Users/ericc/Dropbox/eric chuu research/GGM/Wishart_InvA_RNG.R")
source("C:/Users/ericc/Dropbox/eric chuu research/GGM/HIWsim.R")

#####################################################################
############# HELPER FUNCTIONS FOR CALCULATING LOG MARG LIK #########
#####################################################################

# log multivariate gamma function Gamma_p(a)
logmultigamma = function(p, a){
    f = 0.25*p*(p-1)*log(pi)
    for(i in 1:p){ f = f + lgamma(a+0.5-0.5*i) }
    return(f)
}

logfrac = function(dim, b, D){
    temp = b+dim-1
    logfrac = 0.5*temp*log(det(D)) - 0.5*temp*dim*log(2) - 
        logmultigamma(dim, 0.5*temp)
    return(logfrac)
}

# log normalizing constant for HIW
logHIWnorm = function(G, b, D){
    junct = makedecompgraph(G)
    cliques = junct$C; separators = junct$S
    nc = length(cliques); ns = length(separators)
    
    Cnorm = 0
    for(i in 1:nc){
        ind = cliques[[i]]
        Cnorm = Cnorm + logfrac(length(ind), b, D[ind, ind, drop = FALSE])
    }
    
    Snorm = 0
    if(ns>1){
        for(i in 2:ns){
            ind = separators[[i]]
            Snorm = Snorm + logfrac(length(ind), b, D[ind, ind, drop = FALSE])
        }
    }
    
    logHIWnorm = Cnorm - Snorm
    return(logHIWnorm)
}

# log marginal likelihood log(f(Y|G))
logmarginal = function(Y, G, b, D){
    n = nrow(Y); p = ncol(Y); S = t(Y)%*%Y
    logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) - 
        logHIWnorm(G, b+n, D+S)
    return(logmarginal)
}


# inverse Wishart density
logSingle = function(dim, b, D, Sigma){
    temp = b + dim - 1
    logSingle = 0.5 * temp * log(det(D)) - 0.5 * temp * dim * log(2) - 
        logmultigamma(dim, 0.5 * temp)
    logSingle = logSingle - 0.5 * (b + 2 * dim) * log(det(Sigma)) -
        0.5 * sum(diag( solve(Sigma) %*% D ))
    
    return(logSingle)
}

# HIW density
logHIW = function(G, b, D, Sigma){
    junct = makedecompgraph(G)
    cliques = junct$C; separators = junct$S
    nc = length(cliques); ns = length(separators)
    
    Cnorm = 0
    for(i in 1:nc){
        ind = cliques[[i]]
        Cnorm = Cnorm + logSingle(length(ind), b, 
                                  D[ind, ind, drop = FALSE], 
                                  Sigma[ind, ind, drop = FALSE])
    }
    
    Snorm = 0
    if(ns>1){
        for(i in 2:ns){
            ind = separators[[i]]
            Snorm = Snorm + logSingle(length(ind), b, 
                                      D[ind, ind, drop = FALSE], 
                                      Sigma[ind, ind, drop = FALSE])
        }
    }
    
    logHIW = Cnorm - Snorm
    return(logHIW)
}


