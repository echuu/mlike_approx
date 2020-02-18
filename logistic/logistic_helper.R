

library(numDeriv)


# = -log p(y | u) - log p(u)
psi = function(u, prior) {
    
    y   = prior$y
    X   = prior$X
    D   = prior$D
    tau = prior$tau
    
    Xbeta = X %*% u
    
    # Xbeta = X %*% beta
    
    # y[n] * X[n,], each row of X * corresponding element in y
    # beta %*% colSums(y * X) %>% c # first term
    
    loglik = (y * Xbeta - log(1 + exp(Xbeta))) %>% sum
    logprior = D / 2 * log(tau) - D / 2 * log(2 * pi) - 
        tau / 2 * sum(u * u)
    
    out = - loglik - logprior
    
    return(out)
    
} # end psi() function


# analytic form of the gradient
lambda = function(u, prior) {
    
    
    y   = prior$y
    X   = prior$X
    tau = prior$tau
    
    grad_beta = - colSums(((y - c(inv.logit(X %*% u)))) * X) + tau * u
    
    return(grad_beta %>% c)
}


# compute gradient numerically, compare with analytic gradient
# grad_psi = function(u, prior) {
#     grad(psi, u, prior = prior)
# }
#

# x_test = rnorm(1000)
# 
# library(microbenchmark)
# microbenchmark(
#                "numer" = {
#                    grad_psi(u_samps[1,], prior = prior)
#                },
#                "analytic" = {
#                    lambda(u_samps[1,], prior = prior)
#                })


# f = function(u, prior) {
#     # sum(log(1 + exp(prior$X %*% u)))
#     y   = prior$y
#     X   = prior$X
#     D   = prior$D
#     tau = prior$tau
#     Xbeta = X %*% u
#     # -(y * Xbeta) %>% sum()
#     # -(y * Xbeta + log(1 + exp(Xbeta))) %>% sum
#     (y * Xbeta + log(1 + exp(Xbeta))) %>% sum
# }
# 
# grad_logit = function(u, prior) {
#     grad(f, u, prior = prior)    
# }













