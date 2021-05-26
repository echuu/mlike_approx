


# extract only the free elements from Psi in vector form
extractFree = function(u, free_ind) {
    u[free_ind]
}



# compute Psi = Phi * P^(-1)
computePsi = function(phi_vec, P_inv) {
    p = nrow(P_inv)
    Phi = matrix(phi_vec, p, p, byrow = FALSE)
    Phi %*% P_inv
}



## compute the gwish density
# test = function(u, obj) {
#
#     ## first reconstruct entire Psi matrix
#     Psi_copy = matrix(0, p, p)
#     Phi = matrix(obj$z0, p, p)
#     Psi = matrix(obj$z1, p, p)
#
#     i = 1
#     for (c in 1:p) {
#         for (r in 1:c) {
#             if (G_5[r, c] > 0) {     # if free element, copy contents of u into Psi
#                 Psi_copy[r,c] = u[i]
#                 i = i + 1
#             } else  {         # non-free element
#                 Psi_copy[r,c] = Psi[r,c]
#                 # - sum(Psi_copy[r, r:(s-1)] * P[r:(s-1), s] / P[s, s]) + Phi[r,s] / P[s, s]
#             }
#         }
#     }
#     # Psi
#     # Psi_copy = Psi
#
#     # print(Psi_copy)
#
#     # constant term
#     t0 = sum(0.5 * (b + nu_i) * log(2) + 0.5 * nu_i * log(2 * pi) +
#                  lgamma(0.5 * (b + nu_i))) +
#         sum(0.5 * (b + b_i - 1) * log(diag(P)^2))
#
#     ## compute term that sums over non-free terms
#     NON_FREE = !edgeInd
#     UD_FREE  = (UPPER_DIAG & G_5)
#     diag(UD_FREE) = FALSE
#
#     # non-free terms
#     psi_non_free = Psi_copy[UPPER_DIAG][NON_FREE]
#     t1 = - 0.5 * sum(psi_non_free^2)
#
#     # product over diagonal terms
#     t2 = sum(-lgamma(0.5 * (b + nu_i)) +
#                  (0.5 * (b + nu_i) - 1) *
#                  log(0.5 * diag(Psi_copy)^2) - 0.5 * diag(Psi_copy)^2)
#
#     # product over off-diagonal terms, free terms
#     psi_free_ud = Psi_copy[UD_FREE]
#     t3 = sum(-log(2 * pi) - 0.5 * psi_free_ud^2)
#
#     t_sum = t0 + t1 + t2 + t3
#     return(-t_sum)
#
# }


## sampleGW function: returns the following:
## (1) Psi: free parameters stored as a row vector; this is fed into hybrid algo
## (2) Phi: entire (p x p) matrix where Phi'Phi = Omega, Omega ~ GW(b, V)
samplegw = function(J, G, b, N, V, S, P, param_ind) {

    P_inv = solve(P)

    ## sample from gw distribution: stored as J - (D x D) matrices
    Omega_post = rgwish(J, G_5, b + N, V + S)

    ## Compute Phi (upper triangular), stored column-wise, Phi'Phi = Omega
    Phi = apply(Omega_post, 3, chol) # (p^2) x J

    ## Compute Psi
    Psi = apply(Phi, 2, computePsi, P_inv = P_inv)

    ## Extract free elements of Psi
    Psi_free = apply(Psi, 2, extractFree, free_ind = param_ind)

    out = list(Phi = t(Phi),
               Psi = t(Psi),
               Psi_free = t(Psi_free))

}


sampleParams = function(J, G, b, N, v, S) {
    # array, where each element is a (p x p) precision matrix
    Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
    post_samps = t(apply(Omega_post, 3, process_samps, edgeIndex = edgeIndex))
}


process_samps = function(Omega, edgeIndex){
    Lt = chol(Omega)
    Lt_vec = Lt[upper.tri(Lt, diag = T)]
    Lt_vec[edgeIndex]
}


sampleGW = function(J, edgeIndex, G, b, N, V, S) {

    Omega_post    = vector("list", J) # store posterior samples in matrix form
    # Lt_post       = vector("list", J) # store lower cholesky factor
    # post_samps_0  = matrix(0, J, D_0) # store ENTIRE upper diag in vector form
    # post_samps    = matrix(0, J, D_u) # store NONZERO upper diag in vector form

    Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
    post_samps = t(apply(Omega_post, 3, process_samps, edgeIndex = edgeIndex))
    return(post_samps)
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




# gwish_loglik = function(u, params) {
#
#     N   = params$N
#     D   = params$D
#     S   = params$S
#     b   = params$b       # degree of freedom for G-wishart distribution
#     V   = params$V # scale matrix for G-wishart distribution
#
#     Lt = matrix(0, D, D)              # (D x D) lower triangular matrix
#     Lt[upper.tri(Lt, diag = T)] = u   # populate lower triangular terms
#
#     # recall: Sigma^(-1) = LL'
#     loglik = - 0.5 * N * D * log(2 * pi) + N * log_det(Lt) -
#         0.5 * matrixcalc::matrix.trace(t(Lt) %*% Lt %*% S)
#
#     return(loglik)
# }
#
# gwish_logprior = function(u, params) {
#
#     D   = params$D
#     b   = params$b # degree of freedom for G-wishart distribution
#     V   = params$V # scale matrix for G-wishart distribution
#     nu  = params$nu
#
#     Lt = matrix(0, D, D)              # (D x D) upper triangular matrix
#     Lt[upper.tri(Lt, diag = T)] = u   # populate upper triangular terms
#
#     x = (b - 2) * sum(log(diag(Lt))) -
#         0.5 * matrixcalc::matrix.trace(t(Lt) %*% Lt %*% V) +
#         D * log(2) + sum((nu + 1) * log(diag(Lt)))
#
#     sum((b + nu - 1) * log(diag(Lt))) -
#         0.5 * matrixcalc::matrix.trace(t(Lt) %*% Lt %*% V) +
#         D * log(2)
#
# }


## psi() function  -------------------------------------------------------------
# slow_psi = function(u, params) {
#
#     loglik = gwish_loglik(u, params)
#     logprior = gwish_logprior(u, params)
#
#     return(- logprior - loglik)
#
# } # end of psi() function ------------------------------------------------------
#




# grad = function(u, params) {
#
#     D   = params$D
#     D_0 = params$D_0
#     N   = params$N
#     S   = params$S
#     V   = params$V
#     nu  = params$nu
#     xi  = params$xi
#
#     Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
#     Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
#     Lt[upper.tri(Lt, diag = T)] = u
#
#     grad_mat = - diag((xi + N) / diag(Lt)) + Lt %*% (S + V)
#     grad_vec = grad_mat[upper.tri(grad_mat, diag = T)]
#     # grad_vec[params$edgeInd] # consider only terms that have edge in the graph
#
#     return(grad_vec)
# }
#
# # all.equal(grad(u, params),
# #           pracma::grad(slow_psi, u, params = params))
# #
# #
# # hess(u, params)
# # pracma::hessian(slow_psi, u, params = params)
# # all.equal(hess(u, params),
# #           pracma::hessian(slow_psi, u, params = params))
#
# hess = function(u, params) {
#
#     D     = params$D
#     D_0   = params$D_0
#     S     = params$S
#     V     = params$V
#     nu    = params$nu
#     xi    = params$xi
#     t_ind = params$t_ind
#
#     Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
#     Lt[upper.tri(Lt, diag = T)] = u
#
#     H = matrix(NA, D_0, D_0)
#
#     for (r in 1:D_0) {
#
#         i = t_ind[r, 1] # row of 1st partial
#         j = t_ind[r, 2] # col of 1st partial
#
#         c = r
#         while (c <= D_0) {
#
#             k = t_ind[c, 1] # row of 2nd partial
#             l = t_ind[c, 2] # col of 2nd partial
#
#             if (i != k) {
#                 H[r, c] = H[c, r] = 0
#             } else if (i == j && k == i && l > j) {
#                 H[r, c] = H[c, r] = -S[l, j] - V[l,j]
#             } else if (i == j && j == k && k == l) {
#                 H[r, c] = H[c, r] = -1/Lt[i,i]^2 * (N + xi[i]) - S[i,i] - V[i,i]
#             } else if (i != j && k == i && l == j) {
#                 H[r, c] = H[c, r] = -S[l, j] - V[l,j]
#             } else if (i != j && k == i && l > j) {
#                 H[r, c] = H[c, r] = -S[l, j] - V[l,j]
#             }
#             c = c + 1
#         }
#
#     }
#     return(-H)
# }























