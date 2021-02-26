
## global constants, prior parameters ------------------------------------------
DataObs = y
DataCovariates = X1
PriorMeanMu = rep(0,ncol(X1))
PriorPrecisionMu = 0.01*diag(1,ncol(X1))

n = length(DataObs)
y = DataObs
z = rep(1,n) - DataObs
X = as.matrix(DataCovariates)
d = ncol(X)
mu0 = PriorMeanMu
tau0 = PriorPrecisionMu
choltau0 = chol(tau0)
log.dettau0 = 2*sum(log(diag(choltau0)))
logPi = log(pi)
log2Pi = log(2*pi)
## -----------------------------------------------------------------------------


## pima logistic regression dataset

# define psi, gradient, hessian

setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions


logprior = function(u, params) {
    -0.5*d*log2Pi + 0.5*log.dettau0 - 0.5*t(u - mu0)%*%tau0%*%(u-mu0)
}

loglike = function(u, params) {
    l = X%*%u
    g = exp(l)
    logLikelihood = t(y)%*%l - sum(log(1+g))
    return(logLikelihood)
}

psi = function(u, params) {
    -loglike(u, params)-logprior(u, params)
}


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

params = NULL
u_df_all = preprocess(data.frame(T), d, params)

u_df = tail(u_df_all, 1000)
u_df %>% head

hml_approx = hml_const(1, d, u_df, 1000, params)
hml_approx$const_vec


## do the updated approximation now
lambda = function(u, params) {
    l = X %*% u
    g = exp(l)
    mu = g/(1+g)
    logPosteriorGradient = t(X)%*%(y - mu) - tau0%*%(u - mu0)
    return(-logPosteriorGradient)	
}

hess = function(u, params) {
    H = matrix(nrow=d,ncol=d)
    l = X%*%u
    g = exp(l)
    mu = g/(1+g)
    zmu = 1 - mu
    mu.zmu = mu*zmu
    
    for(i in 1:d){
        for(j in 1:d){
            H[i,j] = -t(S[i,j,])%*%mu.zmu
        }
    }
    H = H - tau0 # prior
    return(-H)
}


l1_norm = function(u, u_0) {
    sum(abs(u - u_0))
}


u_df = u_df_all[sample(1:nrow(u_df_all), 1000),]
u_df %>% head

# ------------------------------------------------------------------------------

D = d
## (2) fit the regression tree via rpart()
u_rpart = rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support) 

#### extension starts here -------------------------------------------------

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

lambda_k = apply(psi_df[,1:D], 1, lambda, params = params)

for (k in 1:K) {
    
    # print(k)
    
    u_k = unname(unlist(psi_df[k,1:D]))
    
    # pracma::grad(psi, u_k, params = params)
    # hess(u_k, params)
    
    # H_k = pracma::hessian(psi, u_k, params = params)
    H_k = hess(u_k, params)
    H_k_inv = chol2inv(chol(H_k))
    # H_k_inv %>% is.symmetric.matrix()
    # H_k_inv = (H_k_inv + t(H_k_inv)) / 2
    # H_k_inv = solve(H_k)
    
    # lambda_k = pracma::grad(psi, u_k, params = params)
    b_k = H_k %*% u_k - lambda_k[,k]
    m_k = H_k_inv %*% b_k
    
    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
    G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])
    
    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
        psi_df$psi_u[k] + sum(lambda_k[,k] * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
        0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
    
}

log_sum_exp(log_terms)














