


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

LIL = -0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) + 
    lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)


params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, beta = beta,
              d = d, n = n, 
              X = X, y = y)







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

# u_0 = u_df[which(u_df$psi_u == max(u_df$psi_u)),1:D] %>% unname() %>% unlist()

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
    
    H_k = hessian(u_k, params = params)
    H_k_inv = chol2inv(chol(H_k))
    
    # lambda_k = pracma::grad(psi, u_k, params = params)
    lambda_k = grad(u_k, params)
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


hybml(u_df, params, psi, grad, hess)

LIL





# grad_psi = function(theta, params) {
#     dist = theta[1:d] - beta0
# 
#     tau = theta[d+1]
# 
#     w = -tau*(M%*%dist)
# 
#     z = 0.5*(n+d)/tau -
#         0.5*t(dist)%*%M%*%dist -
#         0.5*c0 + (0.5*alpha-1)/tau - 0.5*delta
# 
#     return(-c(w,z))
# }
# 
# hess = function(theta, params) {
#     dist = theta[1:d] - beta0
#     tau = theta[d+1]
#     H = matrix(nrow=d+1,ncol=d+1)
#     H[1:d,1:d] = -tau*M
#     z = -M%*%dist
#     H[1:d,d+1] = z
#     H[d+1,1:d] = t(z)
#     H[d+1,d+1] = -0.5*(n+d)/(tau^2) - (0.5*alpha-1)/(tau^2)
#     return(-H)
# }




### testing C++ code
logprior(u1, params)
old_logprior(u1)

loglike(u1, params)
old_loglike(u1, params)

psi(u1, params)
old_psi(u1, params)

library(microbenchmark)
microbenchmark(f1 = logprior(u1, params),
               f2 = old_logprior(u1))

microbenchmark(f1 = loglike(u1, params),
               f2 = old_loglike(u1))

test_samps = u_samps[1:1000,]
test_samps %>% head

test_df_slow = slow_preprocess(test_samps, D)
test_df_slow %>% head

test_df_cpp = preprocess(test_samps, D, params)
test_df_cpp %>% head

microbenchmark(f1 = preprocess(test_samps, D, params),
               f2 = slow_preprocess(test_samps, D))



