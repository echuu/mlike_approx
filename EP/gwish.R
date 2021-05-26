

source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")
library(BDgraph)
# Rcpp::sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
# source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")
# source("C:/Users/ericc/mlike_approx/covariance/HIW_helper.R")

options(scipen = 999)
options(dplyr.summarise.inform = FALSE)

b = 3

d = 5 # number vertices
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), d, d)

edgeInd  = G_5[upper.tri(G_5, diag = TRUE)] %>% as.logical
upperInd = G_5[upper.tri(G_5)] %>% as.logical
D_u = sum(edgeInd)

N = 200
V = diag(1, d)

# Generate data Y based on Sigma_G
set.seed(1)
Omega_G = rgwish(1, G_5, b, V) # generate the true precision matrix

## compute some values used in the log posterior formula
nu_i = numeric(d)
for (i in 1:d) {
  ne = which(G_5[i,] > 0)
  nu_i[i] = length(ne[ne > i])
}

k_i = (1:d) - 1
b_i = nu_i + k_i + 1
b_i

# ------------------------------------------------------------------------------

V   = diag(1, d)


Sigma_G = solve(Omega_G)       # invert for covariance matrix
Y = matrix(0, N, d)            # generate data
for (i in 1:N) {
  Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(d, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

- 0.5 * d * N * log(2 * pi) +
  gnorm(G_5, b + N, V + S, iter = 1000) - gnorm(G_5, b, V, iter = 1000)

nu_i = numeric(d)
for (i in 1:d) {
  ne = which(G_5[i,] > 0)
  nu_i[i] = length(ne[ne > i])
}

xi = b + nu_i - 1

D = d * (d + 1) / 2 # dimension of the parameter (# of diag + upper diag)



# construct D_0 x 2 matrix, where each row has the row, column of each of the
# free parameters **in order** (order matters because this is how they are
# populated in both the gradient vector and the hessian matrix)
index_mat = matrix(0, d, d)
index_mat[upper.tri(index_mat, diag = T)] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

params = list(N = N, D = d, D_0 = D, S = S, b = b, V = V,
              G = G_5, nu = nu_i, xi = xi,
              t_ind = t_ind)

J = 5000
post_gW = sampleGW(J, D, G_5, b, N, V, S) %>% data.frame()


u_df = preprocess(post_gW, D, params)     # J x (D_u + 1)
u_df %>% head()

u = u_df[1,1:D] %>% unlist %>% unname
Lt = matrix(0, d, d)              # (D x D) lower triangular matrix
Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
Lt

## (2) fit the regression tree via rpart()
u_rpart = rpart::rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)

#### extension starts here -------------------------------------------------

### (1) find global mean
# u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean
# u_star = u_0

# grad = function(u, params) { pracma::grad(slow_psi, u, params = params) }
# hess   = function(u, params) { pracma::hessian(slow_psi, u, params = params) }
psi = slow_psi
u_star = globalMode(u_df)
# cbind(u_0, globalMode(u_df))


### (2) find point in each partition closest to global mean (for now)
# u_k for each partition
u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_star)
# cost = apply(u_df_part[,1:D], 1, l2_norm, u_0 = u_star)
u_df_part = u_df_part %>% dplyr::mutate(cost = cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>%
  group_by(leaf_id) %>% filter(cost == min(cost)) %>%
  data.frame

## pick each point whose posterior probability is highest ----------------------
### (1) find global mean
# MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
# u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
# u_star = u_0
#
# psi_df = u_df_part %>% dplyr::group_by(leaf_id) %>%
#   dplyr::filter(psi_u == max(psi_u)) %>%
#   data.frame
# psi_df %>% head
# ------------------------------------------------------------------------------

bounds = u_partition %>% arrange(leaf_id) %>%
  dplyr::select(-c("psi_hat", "leaf_id"))
psi_df = psi_df %>% arrange(leaf_id)
psi_df

K = nrow(bounds)
log_terms = numeric(K) # store terms so that we can use log-sum-exp()
G_k = rep(NA, K)       # store terms coming from gaussian integral


k = 1
source("C:/Users/ericc/rcpp-epmgp/util.R")
for (k in 1:nrow(bounds)) {

  u_k = unname(unlist(psi_df[k,1:D]))

  # H_k = pracma::hessian(slow_psi, u_k, params = params)
  H_k = hess(u_k, params)
  H_k_inv = chol2inv(chol(H_k))

  # lambda_k = pracma::grad(slow_psi, u_k, params = params)
  lambda_k = grad(u_k, params = params)
  b_k = H_k %*% u_k - lambda_k
  m_k = H_k_inv %*% b_k

  lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
  ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

  # source("C:/Users/ericc/rcpp-epmgp/util.R")
  G_k[k] = hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub, EPS_CONVERGE = 1e-5)$logZ
  # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
  # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

  log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
    psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k +
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]

  # k = k + 1

}

log_terms
log_sum_exp(log_terms)

- 0.5 * d * N * log(2 * pi) + BDgraph::gnorm(G_5, b + N, V + S, iter = 1e4) -
  BDgraph::gnorm(G_5, b, V, iter = 1e4)


log_density = function(u, data) {
  -slow_psi(u, data)
}
u_samp = as.matrix(post_gW)
colnames(u_samp) = names(u_df)[1:D]
lb = rep(-Inf, D)
ub = rep(Inf, D)
names(lb) <- names(ub) <- colnames(u_samp)
bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params,
                                               lb = lb, ub = ub,
                                               silent = TRUE)

bridge_result$logml




while (k <= nrow(bounds)) {

  u_k = unname(unlist(psi_df[k,1:D]))

  # H_k = pracma::hessian(slow_psi, u_k, params = params)
  H_k = hess(u_k, params)
  H_k_inv = chol2inv(chol(H_k))

  # lambda_k = pracma::grad(slow_psi, u_k, params = params)
  lambda_k = grad(u_k, params = params)
  b_k = H_k %*% u_k - lambda_k
  m_k = H_k_inv %*% b_k

  lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
  ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

  G_k[k] = hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub, EPS_CONVERGE = 1e-5)$logZ
  # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
  # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

  log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
    psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k +
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]

  k = k + 1

}

G_k
log_terms


log_sum_exp(log_terms)
- 0.5 * d * N * log(2 * pi) +
  gnorm(G_5, delta + N, V + S, iter = 1e4) - gnorm(G_5, delta, V, iter = 1e4)


#### testing gradient and hessian
u = u_df[1,1:D] %>% unlist %>% unname











