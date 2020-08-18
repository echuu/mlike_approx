
# setwd("C:/Users/ericc/Dropbox/eric chuu research/GGM/")

library(BDgraph)
# source("C:/Users/ericc/Dropbox/eric chuu research/GGM/makedecompgraph.R")
# source("C:/Users/ericc/Dropbox/eric chuu research/GGM/Wishart_InvA_RNG.R")
# source("C:/Users/ericc/Dropbox/eric chuu research/GGM/HIWsim.R")



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

# # log marginal likelihood log(f(Y|G))
# logmarginal = function(Y, G, b, D){
#     n = nrow(Y); p = ncol(Y); S = t(Y)%*%Y
#     logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) - 
#         logHIWnorm(G, b+n, D+S)
#     return(logmarginal)
# }

logmarginal = function(Y, G, b, D, S){
    n = nrow(Y); p = ncol(Y); # S = t(Y)%*%Y
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


# Sample the HIW_G(bG,DG) distribution on a graph G with adjacency matrix Adj
HIWsim = function(Adj, bG, DG){
    # check if "Adj" is a matrix object
    if(is.matrix(Adj)==FALSE) { stop("Adj must be a matrix object!") }
    # check if "Adj" is a square matrix
    if(dim(Adj)[1]!=dim(Adj)[2]) { stop("Adj must be a square matrix") }
    # check if "Adj" is symmetric
    if(isSymmetric.matrix(Adj)==FALSE) { stop("Adj must be a symmetric matrix") }
    
    # check if "DG" is a matrix object
    if(is.matrix(DG)==FALSE) { stop("DG must be a matrix object!") }
    # check if "DG" is a square matrix
    if(dim(DG)[1]!=dim(DG)[2]) { stop("DG must be a square matrix") }
    # check if "DG" is symmetric
    if(isSymmetric.matrix(DG)==FALSE) { stop("DG must be a symmetric matrix") }
    
    # check if "bG" is greater than 2
    if(bG<=2) { stop("bG must be greater than 2") }
    # check if "Adj" and "DG" are the same size
    if(nrow(Adj)!=nrow(DG)) { stop("Adj and DG must have the same dimension") }
    rMNorm = function(m, V){ return(m+t(chol(V))%*%rnorm(length(m))) }
    
    p = nrow(Adj)
    temp = makedecompgraph(Adj)
    Cliques = temp$C
    Separators = temp$S
    numberofcliques = length(Cliques)
    
    ############################################################
    # Creat some working arrays that are computed only once
    C1 = solve(DG[Cliques[[1]], Cliques[[1]]]/bG)
    c1 = Cliques[[1]]
    UN = c1
    DSi = DRS = mU = list()
    
    for(i in 2:numberofcliques){
        sid = Separators[[i]]
        DSi[[i]] = solve(DG[sid, sid])
        cid = Cliques[[i]]
        dif = sort(setdiff(cid, UN))
        UN = sort(union(cid, UN)) # no need to sort, just playing safe
        sizedif = length(dif)
        DRS[[i]] = DG[dif, dif] - DG[dif, sid] %*% DSi[[i]] %*% DG[sid, dif]
        DRS[[i]] = ( DRS[[i]] + t(DRS[[i]]) )/2
        mU[[i]] = DG[dif, sid] %*% DSi[[i]]
    }
    
    ############################################################
    # MC Sampling
    UN = c1
    Sigmaj = matrix(0, p, p)
    # sample variance mx on first component
    Sigmaj[c1, c1] = solve(Wishart_InvA_RNG( bG+length(Cliques[[1]])-1, DG[Cliques[[1]], Cliques[[1]]] ))
    
    for(i in 2:numberofcliques){ # visit components and separators in turn
        dif = sort(setdiff(Cliques[[i]], UN))
        UN = sort(union(Cliques[[i]], UN)) # probably no need to sort, just playing safe
        sizedif = length(dif)
        sid = Separators[[i]]
        SigRS = solve(Wishart_InvA_RNG( bG+length(Cliques[[i]])-1, DRS[[i]] ))
        Ui = rMNorm( as.vector(t(mU[[i]])), kronecker(SigRS, DSi[[i]]))
        Sigmaj[dif, sid] = t(matrix(Ui, ncol = sizedif)) %*% Sigmaj[sid, sid]
        Sigmaj[sid, dif] = t(Sigmaj[dif, sid])
        Sigmaj[dif, dif] = SigRS + Sigmaj[dif, sid] %*% solve(Sigmaj[sid, sid]) %*% Sigmaj[sid, dif]
    }
    
    # Next, completion operation for sampled variance matrix
    H = c1
    for(i in 2:numberofcliques){
        dif = sort(setdiff(Cliques[[i]], H))
        sid = Separators[[i]]
        h = sort(setdiff(H, sid))
        Sigmaj[dif, h] = Sigmaj[dif, sid] %*% solve(Sigmaj[sid, sid]) %*% Sigmaj[sid, h]
        Sigmaj[h, dif] = t(Sigmaj[dif, h])
        H = sort(union(H, Cliques[[i]])) # probably no need to sort, just playing safe
    }
    Sigma = Sigmaj
    
    # Next, computing the corresponding sampled precision matrix
    Caux = Saux = array(0, c(p, p, numberofcliques))
    cid = Cliques[[1]]
    Caux[cid, cid, 1] = solve(Sigmaj[cid, cid])
    for(i in 2:numberofcliques){
        cid = Cliques[[i]]
        Caux[cid, cid, i] = solve(Sigmaj[cid, cid])
        sid = Separators[[i]]
        Saux[sid, sid, i] = solve(Sigmaj[sid, sid])
    }
    Omega = rowSums(Caux, dims = 2) - rowSums(Saux, dims = 2)
    
    return(list(Sigma = Sigma, Omega = Omega))
}


# Input:  an adjacency  matrix A of a decomposable graph G
# Output: cell array G containing the cliques and separators of G   
#         nodeIDs and nodenames are optional inputs  

makedecompgraph = function(Adj){
    # first check if "Adj" is a matrix object
    if(is.matrix(Adj)==FALSE) { stop("the input must be a matrix object!") }
    # check if "Adj" is a square matrix
    if(dim(Adj)[1]!=dim(Adj)[2]) { stop("the input must be a square matrix") }
    # check if "Adj" is symmetric
    if(isSymmetric.matrix(Adj)==FALSE) { stop("the input must be a symmetric matrix") }
    
    p = nrow(Adj)
    Adj[Adj!=0] = 1 # set all non-zero entries of Adj to 1
    Adj = Adj - diag(diag(Adj)) + diag(p) # set all diagonal entries of Adj to be 1
    Order = 1:p
    i = 1
    Adj0 = Adj
    while(i<p){
        nn = apply(Adj[1:i, (i+1):p, drop=F], 2, sum)
        b = which.max(nn)
        Order[c(i+1, b+i)] = Order[c(b+i, i+1)]
        i = i + 1
        Adj = Adj0
        Adj = Adj[Order, Order]
    }
    
    numberofcliques = 1
    Cliques = list(1)
    i = 2
    while(i<=p){
        if( sum(Adj[i,Cliques[[numberofcliques]]])==length(Cliques[[numberofcliques]]) ){
            Cliques[[numberofcliques]] = c(Cliques[[numberofcliques]], i)
        }
        else{
            numberofcliques = numberofcliques + 1
            Cliques[[numberofcliques]] = union(i, which(Adj[i, 1:i]==1))
        }
        i = i + 1
    }
    
    for(i in 1:numberofcliques){
        Cliques[[i]] = sort(Order[Cliques[[i]]])
    }
    
    UN = Cliques[[1]]
    Separators = list()
    if(numberofcliques==1){ return(list(C = Cliques, S = Separators)) }
    else{
        for(i in 2:numberofcliques){
            Separators[[i]] = sort(intersect(UN, Cliques[[i]]))
            UN = union(UN, Cliques[[i]])
        }
        return(list(C = Cliques, S = Separators))
    }
}



# Generates 1 draw from a Wishart distribution - allows for singular Wishart too, 
# in cases that the distn parameter S is rank deficient of rank r.

# K is W_p(df,A) with sum-of-squares parameter S = A^{-1} and d.o.f. df
# Dimension p is implicit

# Usual nonsingular case: r=p<df; now allow for r<p  with r integral
# Note that  E(K)=df.S^{-1}  in this notation 

# Noticing pdf is p(K) = cons. |K|^((df-p-1)/2) exp(-trace(K S)/2) 
# in case usual modification 

# Returns matrix W of dimension p by p of rank r=rank(S)

# EXAMPLE: reference posterior in a normal model N(0,Sigma) with precision 
# matrix Omega = Sigma^{-1}
# Random sample of size n has sample var matrix V=S/n with S=\sum_i x_ix_i'
# Ref posterior for precision mx Omega is W_p(n,A) with A=S^{-1}
# e.g., K = wishart_InvA_rnd(n,n*V); draw 1000 samples

# Useful for looking at posteriors for correlations and eigenvalues
# and also for posterior on graph - looking for elements of Omega near 0

# K = Wishart_InvA_RNG(df, S) <=> K ~ Wishart_p(df, S^{-1})
# E[K] = df * S^{-1}
Wishart_InvA_RNG = function(df, S){
    # first check if "S" is a matrix object
    if(is.matrix(S)==FALSE) { stop("the input must be a matrix object!") }
    # check if "S" is a square matrix
    if(dim(S)[1]!=dim(S)[2]) { stop("the input must be a square matrix") }
    # check if "S" is symmetric
    if(isSymmetric.matrix(S)==FALSE) { stop("the input must be a symmetric matrix") }
    
    p = nrow(S)
    temp = svd(S)
    P = temp$u
    D = diag(1/sqrt(temp$d), nrow = length(temp$d))
    i = which(temp$d>max(temp$d)*1e-9)
    r = length(i)
    P = P[,i] %*% D[,i]
    U = matrix(0, nrow = p, ncol = r)
    for(j in 1:r){
        U[j, j:r] = c(sqrt( rgamma(1, shape = (df-j+1)/2, scale = 2) ), rnorm(r-j))
    }
    U = U %*% t(P)
    K = t(U) %*% U
    return(K)
}







