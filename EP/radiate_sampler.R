
source("ModelEvidence.R")

n.its = 505000
burn.in = 101000
E1 = evidence.obj(y,X,mu0,Lambda0,3,2*300^2)

# define fix


test = E1$log.harmonic.mean.estimator.evidence(Its = n.its, BurnIn = burn.in)


T = E1$gibbs.sampler(Its,BurnIn,fix)

logLikelihood = apply(t(T),2,function(x) E1$log.likelihood(x))

# try to get stability
a = min(logLikelihood)

harmonic.mean.estimator = log(n) + a - log(sum(exp(-logLikelihood + a)))



gibbs_radiata = function(Its, BurnIn, fix, initial = NULL,
                         return.log.posterior = FALSE,
                         return.log.likelihood = FALSE) {
    
    # do site by site updates for fair comparison between methods
    T = matrix(nrow = Its - BurnIn,ncol = d+1)
    
    #inialize from prior
    if (is.null(initial)) {
        tau = rgamma(1,shape = alpha/2,rate = delta/2)
        beta = rnorm(d,mean = mu0,sd = sqrt(1/(tau*diag(tau0))))
    } else{
        beta = initial[1:d]
        tau = intial[d+1]
    }
    
    sh = 0.5*(n+d+alpha)
    
    sample.vars = which(fix$vars[1:d] == FALSE)
    
    fix.vars = which(fix$vars[1:d] == TRUE)
    if(length(fix.vars)>0) beta[fix.vars] = fix$values[fix.vars]
    
    if(fix$vars[d+1] == TRUE) tau = fix$values[d+1]
    
    sample.tau = !fix$vars[d+1]
    
    for(ItNum in 1:Its){
        
        #visit each parameter in turn
        
        for(j in sample.vars){
            
            w = M[j,]%*%(beta-beta0) - M[j,j]*(beta[j]-beta0[j])
            
            mu = beta0[j] - w/M[j,j]
            
            sig = sqrt(1/(tau*M[j,j]))
            
            beta[j] = rnorm(1,mean=mu,sd=sig)
            
            
        }
        
        rt = 0.5*( t(beta-beta0)%*%M%*%(beta-beta0) + c0 + delta )
        
        if(sample.tau) tau = rgamma(1,shape=sh,rate = rt)
        
        if(ItNum > BurnIn){
            T[ItNum-BurnIn,] = c(beta,tau)
        }
    }		
    
    return(T)
}

