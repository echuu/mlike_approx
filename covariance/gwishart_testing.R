
### gwishart test file

## test equivalence between the numerical gradient, closed form gradient in R,
## and the c++ implementation of the gradient

library(BDgraph)
library(Rcpp)

setwd("C:/Users/ericc/mlike_approx/algo")

source("setup.R")           # setup global environment, load in algo functions
source("C:/Users/ericc/mlike_approx/covariance/gWish_helper.R")


# compute normalizing constant for G_5

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

N = 500
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
N = 500
V = diag(1, D)
# S = t(Y) %*% Y
delta = b


J = 2000
N = 0
S = 0

xi = b + nu_i - 1

params = list(N = N, D = D, D_0 = D_0, S = S, b = b, V = V, 
              G = G_5, nu = nu_i, xi = xi)

post_gW = sampleGW(J, D_0, G_5, b, N, V, S) %>% data.frame()

sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")

u_df = preprocess(post_gW, D_0, params)     # J x (D_u + 1) 
u_df %>% head()

# ------------------------------------------------------------------------------
## begin testing the gradient

u = u_df[1,1:D_0] %>% unlist %>% unname
Lt = matrix(0, p, p)              # (D x D) upper triangular matrix
Lt[upper.tri(Lt, diag = T)] = u   # populate upper triangular terms
Lt

# numerical gradient
test_grad = pracma::grad(slow_psi, u, params = params)
test_grad
grad_mat = matrix(0, p, p)              # (D x D) upper triangular matrix
grad_mat[upper.tri(grad_mat, diag = T)] = test_grad   # populate upper triangular terms
grad_mat

# closed form implementation of gradient in R
closed_grad = - diag(xi/diag(Lt)) + Lt %*% V
closed_grad[upper.tri(closed_grad, diag = T)]

# check equality between numerical and closed form gradient
all.equal(grad_mat, closed_grad)

# load the C++ implementation of the gradient
sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
grad_cpp = grad_gwish(u, params)

# check equality between closed from gradient (R) and C++ implementation
all.equal(closed_grad[upper.tri(closed_grad, diag = T)], grad_cpp)


# ------------------------------------------------------------------------------
## begin testing the hessian

H = pracma::hessian(slow_psi, u, params = params)
diag(H)


#### final version of the hessian ----------------------------------------------
h = function() {
    G = matrix(NA, D_0, D_0) # populate G with: d t_{i,j} d t_{k, l}
    for (r in 1:D_0) {
        i = t_ind[r,1]
        j = t_ind[r,2]
        c = r
        while (c <= D_0) {
            k = t_ind[c,1] # row of 2nd order partial
            l = t_ind[c,2] # col of 2nd order partial
            # partial derivative cases
            MISMATCH_PARTIAL = (i != k)
            DIAGONAL_PARTIAL = (i == j && j == k && k == l)
            MIXED_PARTIAL    = (i != j && k == i && l >= j)
            
            if (MISMATCH_PARTIAL) {
                G[r,c] = 0
            } else if (DIAGONAL_PARTIAL) { 
                G[r,c] = xi[i] / Lt[i,i]^2 + V[i,i]
            } else if (MIXED_PARTIAL) {
                G[r,c] = V[l,j]
            }
            c = c + 1
        }
        
    }
    G[is.na(G)] = 0
    G = 0.5 * (G + t(G))
    G
}
#### final version of the hessian ----------------------------------------------

h()
H = pracma::hessian(slow_psi, u, params = params)
H
all.equal(h(), H)


sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
hess_gwish(u, params)

all.equal(h(), H)
all.equal(hess_gwish(u, params), H)

diag(h())
diag(hess_gwish(u, params))
diag(H)


library(microbenchmark)
microbenchmark(numer = pracma::hessian(slow_psi, u, params = params), 
               r = h(),
               cpp = hess_gwish(u, params))

#### test vec2chol() function --------------------------------------------------

sourceCpp("C:/Users/ericc/mlike_approx/speedup/gwish.cpp")
v2c = function(u, D) {
    Lt = matrix(0, D, D)              # (D x D) upper triangular matrix
    Lt[upper.tri(Lt, diag = T)] = u   # populate upper triangular terms
    Lt
}

vec2chol(u, D)
v2c(u, D)


microbenchmark(r = v2c(u, D),
               cpp = vec2chol(u, D))


#### test chol2vec() function --------------------------------------------------
all.equal(Lt[upper.tri(Lt, diag = T)], chol2vec(Lt, D))
microbenchmark(r = Lt[upper.tri(Lt, diag = T)],
               cpp = chol2vec(Lt, D))










