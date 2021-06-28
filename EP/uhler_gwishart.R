
## gwishart using uhler way of computing the normalizing constant


setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           # setup global environment, load in algo functions


source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")
library(BDgraph)
# library(Rcpp)

# Eq. (2.4) in Uhler paper (p.12)



d = 5
delta = 100
# |E| = 7. Applying formula at the end of pg.2
# of the Uhler paper.
#log normalizing  constant of the G-WIshart is
(LIL = log(2^(0.5*d*delta + 7)) + I_G(0.5*(delta-2)))


G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), d, d)

V = diag(1, d)
gnorm(G_5, delta, V, 100)
J = 1000
N = 0
S = matrix(0, d, d)
post_gW = sampleGW(J, D_0, G_5, delta, N, V, S) %>% data.frame()


nu_i = numeric(d)
for (i in 1:d) {
    ne = which(G_5[i,] > 0)
    nu_i[i] = length(ne[ne > i])
}
nu_i
xi = delta + nu_i - 1
xi
# construct D_0 x 2 matrix, where each row has the row, column of each of the
# free parameters **in order** (order matters because this is how they are
# populated in both the gradient vector and the hessian matrix)
p = d
D = d * (d + 1) / 2
index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

params = list(N = N, D = d, D_0 = D, S = S, b = delta, V = V,
              G = G_5, nu = nu_i, xi = xi,
              t_ind = t_ind)


## create separate helper file for uhler example

u_df = preprocess(post_gW, D, params)     # J x (D_u + 1)
u_df %>% head()


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
gnorm(G_5, delta, V, 1000)


log_density = function(u, data) {
    -psi(u, data)
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
                                               silent = TRUE,
                                               method = 'warp3')

bridge_result$logml






# compute normalizing constant for G_5 -----------------------------------------

# Eq. (2.4) in Uhler paper (p.12)
I_G = function(delta) {
    7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
        lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
        lgamma(delta + 5/2)
}


delta = 3
I_G(delta)

p = 5 # number vertices
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)


N = 10
V = diag(1, p)
delta_post = delta + N

I_G(delta_post)
gnorm(G_5, delta_post, V, 100)

is.symmetric.matrix(G_5)
D = p
nu_i = numeric(D)
for (i in 1:D) {
    ne = which(G_5[i,] > 0)
    nu_i[i] = length(ne[ne > i])
}
nu_i

# target quantity is Eq. (2.2) on p.9 of Uhler
p = 5 # number vertices
D = p
D_0 = D * (D + 1) / 2
b = delta_post
# N = 500
V = diag(1, D)
# S = t(Y) %*% Y
delta = b


J = 10000
N = 0
S = 0

xi = b + nu_i - 1

# construct D_0 x 2 matrix, where each row has the row, column of each of the
# free parameters **in order** (order matters because this is how they are
# populated in both the gradient vector and the hessian matrix)
index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)] = 1:D_0
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

params = list(N = N, D = D, D_0 = D_0, S = S, b = b, V = V,
              G = G_5, nu = nu_i, xi = xi,
              t_ind = t_ind)

post_gW = sampleGW(J, D_0, G_5, b, N, V, S) %>% data.frame()

Rcpp::sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")

source("C:/Users/ericc/mlike_approx/EP/hybml.R")

u_df = preprocess(post_gW, D_0, params)     # J x (D_u + 1)
u_df %>% head()

u = u_df[1,1:D_0] %>% unlist %>% unname
Lt = matrix(0, p, p)              # (D x D) upper triangular matrix
Lt[upper.tri(Lt, diag = T)] = u   # populate upper triangular terms
Lt
#
# hml_approx = hml_const(1, D_0, u_df, J, params)
# hml_approx$const_vec
#
gnorm(G_5, delta_post, V, 100)
hybridml::hybml_const(u_df)$logz


#
# hybml(u_df, params = params, psi = psi, grad = grad_gwish, hess = hess_gwish)
# ------------------------------------------------------------------------------

#### start here


source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")
library(BDgraph)
Rcpp::sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")


d = 5
delta = 3

p = 5 # number vertices
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)


N = 100

V = diag(1, p)
delta_post = delta + N


# Generate data Y based on Sigma_G
Y = matrix(0, N, d)
for (i in 1:N) {
    Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(d, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

- 0.5 * d * N * log(2 * pi) +
    gnorm(G_5, delta + N, V + S, iter = 1000) - gnorm(G_5, delta, V, iter = 1000)

nu_i = numeric(d)
for (i in 1:d) {
    ne = which(G_5[i,] > 0)
    nu_i[i] = length(ne[ne > i])
}
nu_i

# target quantity is Eq. (2.2) on p.9 of Uhler
p = 5 # number vertices
D = p
D_0 = D * (D + 1) / 2
b = delta
V = diag(1, D)
S = t(Y) %*% Y

J = 10000
# N = 0
# S = 0

xi = b + nu_i - 1

# construct D_0 x 2 matrix, where each row has the row, column of each of the
# free parameters **in order** (order matters because this is how they are
# populated in both the gradient vector and the hessian matrix)
index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)] = 1:D_0
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

params = list(N = N, D = D, D_0 = D_0, S = S, b = b, V = V,
              G = G_5, nu = nu_i, xi = xi,
              t_ind = t_ind)

post_gW = sampleGW(J, D_0, G_5, b, N, V, S) %>% data.frame()


u_df = preprocess(post_gW, D_0, params)     # J x (D_u + 1)
u_df %>% head()

D = D_0
# D = D_u
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

bounds = u_partition %>% arrange(leaf_id) %>%
    dplyr::select(-c("psi_hat", "leaf_id"))
psi_df = psi_df %>% arrange(leaf_id)

K = nrow(bounds)
log_terms = numeric(K) # store terms so that we can use log-sum-exp()
G_k = rep(NA, K)       # store terms coming from gaussian integral

# lambda_k = apply(psi_df[,1:D], 1, lambda, params = params)

k = 3

for (k in 1:nrow(bounds)) {

    u_k = unname(unlist(psi_df[k,1:D]))

    H_k = pracma::hessian(slow_psi, u_k, params = params)
    # H_k = hess(u_k, params)
    H_k_inv = chol2inv(chol(H_k))

    lambda_k = pracma::grad(psi, u_k, params = params)
    # lambda_k = grad(u_k, params = params)
    b_k = H_k %*% u_k - lambda_k
    m_k = H_k_inv %*% b_k

    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

    # epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # epmgp::axisepmgp(m_k, H_k_inv, lb, ub)

    # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # G_k[k] = ep_step(lb, ub, m_k, H_k_inv)

    # epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)

    # hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub, EPS_CONVERGE = 1e-5)$logZ

    G_k[k] = hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub, EPS_CONVERGE = 1e-5)$logZ
    # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
        psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k +
        0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]

}

G_k
log_terms


log_sum_exp(log_terms)

- 0.5 * d * N * log(2 * pi) +
    gnorm(G_5, delta + N, V + S, iter = 1e4) - gnorm(G_5, delta, V, iter = 1e4)

# mean(G_k)
# G_k[G_k == 0] = mean(G_k)
#
# log_terms[which(G_k == 0)] = log_terms[which(G_k == 0)] + mean(G_k[G_k!=0])
# log_terms

log_sum_exp(log_terms - G_k)
log_sum_exp(log_terms)
gnorm(G_5, delta, V, 1000)


log_density = function(u, data) {
    -slow_psi(u, data)
}
u_samp = as.matrix(post_gW)
colnames(u_samp) = names(u_df)[1:D_0]
lb = rep(-Inf, D_0)
ub = rep(Inf, D_0)
names(lb) <- names(ub) <- colnames(u_samp)
bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params,
                                               lb = lb, ub = ub,
                                               silent = TRUE,
                                               method = 'warp3')

bridge_result$logml



for (k in 1:K) {

    u_k = unname(unlist(psi_df[k,1:D]))

    H_k = pracma::hessian(psi, u_k, params = params)
    H_k_inv = chol2inv(chol(H_k))

    if (!is.positive.definite(H_k)) {
        print(k)
    }
}


# ------------------------------------------------------------------------------




library(matlabr)
options(matlab.path = "C:/Users/ericc/Dropbox/epmgp/")
have_matlab()

code = c("x = 10",
         "y=20;",
         "z=x+y",
         "a = [1 2 3; 4 5 6; 7 8 10]",
         "save('test.txt', 'x', 'a', 'z', '-ascii')")
res = run_matlab_code(code)


test

H_k_inv

library(R.matlab)
writeMat(con = "C:/Users/ericc/Dropbox/epmgp/hmat.mat", x = as.matrix(H_k_inv))
writeMat(con = "C:/Users/ericc/Dropbox/epmgp/m_k.mat", x = as.matrix(m_k))
writeMat(con = "C:/Users/ericc/Dropbox/epmgp/lb.mat", x = as.matrix(lb))
writeMat(con = "C:/Users/ericc/Dropbox/epmgp/ub.mat", x = as.matrix(ub))








