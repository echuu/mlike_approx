


k = 2

u_k = unname(unlist(psi_df[k,1:D]))

# H_k = pracma::hessian(slow_psi, u_k, params = params) # numerical hessian
H_k = hess_gwish(u_k, params) # closed form computation of hessian -> precision
H_k_inv = chol2inv(chol(H_k)) # inverse precision -> covariance

# lambda_k = pracma::grad(slow_psi, u_k, params = params) # numerical gradient
lambda_k = grad_gwish(u_k, params) # closed form computation of gradient

b_k = H_k %*% u_k - lambda_k
m_k = H_k_inv %*% b_k

lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

# this next line fails
epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)


G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)

log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - 
    psi_df$psi_u[k] + sum(lambda_k[,k] * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k + 
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]



