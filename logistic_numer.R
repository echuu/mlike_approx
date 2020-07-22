vblogit <- function(y, X, offset, eps=1e-2, m0, S0, S0i, xi0, verb=FALSE, maxiter=1000, ...) {
    ### Logistic regression using JJ96 idea.
    ## p(y, w, t) = p(y | w) p(w | t) p(t) 
    ##
    ## Y ~ Bern(logit(Xw + offset))
    ## w  ~ N(m0, S0) iid
    ##
    ##
    cat2 <- if(verb) cat else function(...) NULL
    varnames <- colnames(data.frame(as.matrix(X[1:2,])))
    
    ## Write 
    X <- Matrix(X)
    X <- drop0(X)
    N <- length(y)
    K <- ncol(X)
    #
    #
    # offset
    if(missing('offset')) offset <- 0
    if(length(offset)<N) offset <- rep(offset, N)[1:N]
    #
    #
    # Priors and initial estimates.
    if(missing(S0))S0   <- Diagonal(1e5, n=K)
    if(missing(S0i))S0i <- solve(S0)
    if(missing(m0))m0   <- rep(0, K)
    #' Constants:
    oo2 <- offset^2
    LE_CONST <- as.numeric( -0.5*t(m0)%*%S0i%*%m0 - 0.5*determinant(S0)$mod + sum((y-0.5)*offset) ) 
    Sm0 <- S0i%*%m0
    #' start values for xi:
    if(missing(xi0))xi0   <- rep(4, N)
    if(length(xi0)!=N) xi0 <- rep(xi0, N)[1:N]
    
    est <- list(m=m0, S=S0, Si=S0i, xi=xi0)
    #'
    #
    ## helper functions needed:
    Afunc <- function(x)  -tanh(x/2)/(4*x)
    Cfunc <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
    ###
    ## loop
    le <- -Inf
    le_hist <- le
    loop <- TRUE
    iter <- 0
    #' initials:
    la <- Afunc(xi0)
    Si <- S0i - 2 * t(X*la)%*%X
    S <- solve(Si)
    m <- S%*%( t(X)%*%( (y-0.5) + 2*la*offset ) + Sm0  )
    #
    elboplot <-numeric(maxiter)
    ctr=1
    # Main loop:
    while(loop){
        old <- le
        #' update variational parameters
        M <- S+m%*%t(m)
        #' force symmetric in case of tiny numerical errors
        M <- (M+t(M))/2
        L <- t(chol(M))
        V <- X%*%L
        dR <- rowSums(V^2)
        dO <- 2*offset*(X%*%m)[,1]
        xi2 <- dR + dO + oo2
        
        xi <- sqrt(xi2)
        la <- Afunc(xi)
        # update post covariance
        Si <- S0i - 2 * t(X*la)%*%X
        S <- solve(Si)
        # update post mean
        m <- S%*%( t(X)%*%( (y-0.5) + 2*la*offset ) + Sm0  )
        ## compute the log evidence
        le <-  as.numeric( 0.5*determinant(S)$mod + sum( Cfunc(xi) ) + sum(oo2*la) + 0.5*t(m)%*%Si%*%m + LE_CONST    )
        elboplot[ctr] <-le
        ctr=ctr+1
        # check convergence
        d <- le - old
        if(d < 0) warning("Log-evidence decreasing, try different starting values for xi.")
        loop <- abs(d) > eps & (iter<-iter+1) <= maxiter
        le_hist <- c(le_hist, le)
        cat2("diff:", d, "             \r")
    }
    elboplot=elboplot[1:iter]
    if(iter == maxiter) warning("Maximum iteration limit reached.")
    cat2("\n")
    ## done. Compile:
    est <- list(m=m, S=S, Si=Si, xi=xi, lambda_xi=la)
    # 
    # Marginal evidence
    est$logLik <- le
    est$elboplot<- elboplot
    #
    # Compute max logLik with the Bernoulli model, this should be what glm gives:
    est$logLik_ML <- as.numeric( t(y)%*%(X%*%m+offset) - sum( log( 1 + exp(X%*%m+offset)) ) )
    # 
    # Max loglik with the approximation
    est$logLik_ML2 <- as.numeric(  t(y)%*%(X%*%m + offset)  + t(m)%*%t(X*la)%*%X%*%m - 0.5*sum(X%*%m) + sum(Afunc(xi)) +
                                       2*t(offset*la)%*%X%*%m + t(offset*la)%*%offset - 0.5 * sum(offset)  )
    # 
    # some additional parts, like in glm output
    est$coefficients <- est$m[,1]
    names(est$coefficients) <- varnames
    est$converged <- !(maxiter==iter)
    # additional stuff
    est$logp_hist <- le_hist
    est$parameters <- list(eps=eps, maxiter=maxiter)
    est$priors <- list(m=m0, S=S0)
    est$iterations <- iter
    
    ## return
    
    est
}



library(VGAM)
library(dplyr)
library(Matrix)

set.seed(1)
N = 50
D = 2
X = matrix(rnorm(N * D, 0, sqrt(1/N)), nrow = N)
beta = runif(D, 0.5, 3) # * sign(runif(D, -1,1)))
z = X %*% beta
pr = exp(z) / (1 + exp(z))
y = rbinom(N, 1, pr)

# df = data.frame(y, X)

## prior params
tau = 1
prior = list(y = y, X = X, D = D, N = N, tau = tau)


# run VB to get posterior mean, covariance estimate
fit_vb <- vblogit(y, X, verb=TRUE) 
mu_beta = fit_vb$m # VB posterior mean
Sig_beta = matrix(fit_vb$S, D, D) # VB posterior covariance


# function to numerically integrate
fun_logit = function(b1, b2) {
    
    # Xbeta = X %*% c(b1, b2)
    
    # matrix multiplicatino computed like this because using the above line of
    # code causes complaints from integral2() function
    Xbeta = c(rep(b1, N) * X[,1] + rep(b2, N) * X[,2])
    Xbeta_0 = c(X %*% as.matrix(mu_beta))
    
    # subtract off loglikelihood of beta_hat (VB estimate)
    loglik = (y * (Xbeta - Xbeta_0) - log1pexp(Xbeta) +
                  log1pexp(Xbeta_0)) %>% sum
    
    logprior = D / 2 * log(tau) - D / 2 * log(2 * pi) - 
        tau / 2 * (b1^2 + b2^2)
    
    out = loglik + logprior
    exp(out)
}


library(pracma)

Xbeta_0 = c(X %*% as.matrix(mu_beta))
loglik_0 = (y * Xbeta_0 - log(1 + exp(Xbeta_0))) %>% sum

# numerical integration works for smaller range for limits of integration
result = integral2(fun_logit, -2, 2, -5, 15, reltol = 1e-5)
result$Q
loglik_0 + log(result$Q) 

# increasing the range for the limits of integration makes the integral too small
# log evaluates to 0
result = integral2(fun_logit, -10, 10, -10, 10, reltol = 1e-5)
result$Q
loglik_0 + log(result$Q) 













