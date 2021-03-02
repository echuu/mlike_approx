

setwd("C:/Users/ericc/normalizingconstant/Code Jason Wyse")
source("ModelEvidence.R")

RadiataPine = read.table("RadiataPine.txt",sep = " ",header = TRUE)
y = RadiataPine$y
n = length(y)
X1 = cbind(rep(1,n),RadiataPine$x-mean(RadiataPine$x))


X = X1 # (42 x 2)

mu0 = c(3000,185)
Lambda0 = diag(1,2)
Lambda0[1,1] = 0.06
Lambda0[2,2] = 6.00

a_0 = 3
b_0 = 2*300^2

tX = t(X)
d = ncol(X)
tau0 = Lambda0
alpha = 2 * a_0
delta = 2 * b_0


# n = length(DataObs)
# y = DataObs
# X = as.matrix(DataCovariates)
# tX = t(X)
# d = ncol(X)
# mu0 = PriorMeanMu
# tau0 = PriorSignalMu
# alpha = 2*PriorShapePrecision
# delta = 2*PriorRatePrecision

# convenient constants
logPi = log(pi)
log2Pi = log(2*pi)
XTX = tX%*%X
XTy = tX%*%y
M = XTX + tau0
cholM = chol(M)
log.detM = 2*sum(log(diag(cholM)))
invM = chol2inv(cholM)
choltau0 = chol(tau0)
invtau0 = chol2inv(choltau0)
log.dettau0 = 2*(sum(log(diag(choltau0))))
P = diag(1,n) - X%*%invM%*%tX
beta0 = invM%*%(tX%*%y + tau0%*%mu0)
yTy = t(y)%*%y
c0 = yTy + t(mu0)%*%(tau0%*%mu0) - t(beta0)%*%M%*%beta0
c1 = t(mu0)%*%tau0%*%mu0

-0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) + 
    lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)




params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, beta = beta,
              d = d, n = n, 
              X = X, y = y)

grad = function(u, params) {
    beta = u[1:d] %>% unlist %>% unname
    tau = u[d+1]
    
    g1 = tau * (Lambda0 %*% (beta - mu0) - t(X) %*%(y - X %*% beta))
    g2 = -1/tau * (0.5 * (n + d + alpha) - 1) + 
                       0.5 * (delta + sum((y - X%*%beta)^2) + 
                                  t(beta - mu0) %*% Lambda0 %*% (beta - mu0))
    return(c(g1, g2))
}

u1
grad(u1, NULL)
grad_psi(u1, NULL)

all.equal(grad(u1), grad(u1))


hessian = function(u, params) {
    
    beta = u[1:d] %>% unlist %>% unname
    tau = u[d+1]
    
    h11 = tau * M
    h12 = M %*% beta - XTy - Lambda0 %*% mu0
    h22 = tau^(-2) * (0.5 *(n + d + alpha) - 1)
    
    H = matrix(0, nrow = d + 1, ncol = d + 1)
    H[1:d, 1:d] = h11
    H[d+1, 1:d] = t(h12)
    H[1:d, d+1] = h12
    H[d+1, d+1] = h22
    
    return(H)
    
}

hess(u1)
hessian(u1)

all.equal(hess(u1), hessian(u1))

# ------------------------------------------------------------------------------


# define psi, gradient, hessian
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions

logprior = function(theta, params) {
    dist = theta[1:d] - beta0
    tau = theta[d+1]
    logPrior = -0.5*(d)*log2Pi + 0.5*d*log(tau)+0.5*log.dettau0 - 0.5*tau*t(dist)%*%tau0%*%dist + dgamma(tau,shape = 0.5*alpha,rate = 0.5*delta,log=TRUE)
    
    return(logPrior)
}

loglike = function(theta, params) {
    beta = theta[1:d]
    tau = theta[d+1]
    z = y - X%*%beta
    logLikelihood = -0.5*n*log2Pi + 0.5*n*log(tau) - 0.5*tau*t(z)%*%z
    return(logLikelihood)
}

logpost = function(theta, params) {
    dist = theta[1:d] - beta0
    tau = theta[d+1]
    logPosterior = -0.5*(n+d)*log2Pi + 0.5*(n+d)*log(tau) + 0.5*log.dettau0 -
        0.5*tau*(t(dist)%*%M%*%dist) - 0.5*tau*c0 +0.5*alpha*log(0.5*delta) - 
        lgamma(0.5*alpha) + (0.5*alpha-1)*log(tau) - 0.5*delta*tau
    logPosterior
}


psi = function(u, params) {
    -loglike(u, params)-logprior(u, params)
}

grad_psi = function(theta, params) {
    dist = theta[1:d] - beta0
    
    tau = theta[d+1]
    
    w = -tau*(M%*%dist)
    
    z = 0.5*(n+d)/tau - 
        0.5*t(dist)%*%M%*%dist - 
        0.5*c0 + (0.5*alpha-1)/tau - 0.5*delta
    
    return(-c(w,z))
}

hess = function(theta, params) {
    dist = theta[1:d] - beta0
    tau = theta[d+1]	
    H = matrix(nrow=d+1,ncol=d+1)
    H[1:d,1:d] = -tau*M
    z = -M%*%dist
    H[1:d,d+1] = z
    H[d+1,1:d] = t(z)
    H[d+1,d+1] = -0.5*(n+d)/(tau^2) - (0.5*alpha-1)/(tau^2)
    return(-H)
}


preprocess = function(post_samps, D, params = NULL) {
    
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


n.its = 505000
burn.in = 101000
fix = list(); fix$vars = rep(FALSE,d+1); fix$values = numeric(d+1);
u_samps = gibbs_radiata(Its,BurnIn,fix)
u_samps = u_samps %>% as.data.frame # 200000 x 3
u_samps %>% dim

D = d + 1
u_df_all = preprocess(u_samps, D, params)
row.names(u_df_all) = NULL
u_df_all %>% head
u_df = u_df_all
u_df = u_df_all[sample(1:nrow(u_df_all), 5000),]
u_df %>% head

hml_approx = hml_const(1, D, u_df, 1000, params)
hml_approx$const_vec

# ------------------------------------------------------------------------------
## (2) fit the regression tree via rpart()
u_rpart = rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support) 


### (1) find global mean
u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean

### (2) find point in each partition closest to global mean (for now)
# u_k for each partition
u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>% 
    group_by(leaf_id) %>% filter(l1_cost == min(l1_cost)) %>% 
    data.frame

bounds = u_partition %>% arrange(leaf_id) %>% 
    dplyr::select(-c("psi_hat", "leaf_id")) 
psi_df = psi_df %>% arrange(leaf_id)

K = nrow(bounds)
log_terms = numeric(K) # store terms so that we can use log-sum-exp()
G_k = numeric(K)       # store terms coming from gaussian integral

# lambda_k = apply(psi_df[,1:D], 1, lambda, params = params)

k = 1
for (k in 1:K) {
    
    u_k = unname(unlist(psi_df[k,1:D]))
    
    H_k = hess(u_k, params = params)
    H_k_inv = chol2inv(chol(H_k))
    
    # lambda_k = pracma::grad(psi, u_k, params = params)
    lambda_k = grad_psi(u_k, params)
    b_k = H_k %*% u_k - lambda_k
    m_k = H_k_inv %*% b_k
    
    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
    
    G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    
    # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])
    
    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
        psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
        0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
    
}
log_sum_exp(log_terms)


