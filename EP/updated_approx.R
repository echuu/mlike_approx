



grad = function(u, params) {
    return(params$Q_beta %*% unname(unlist(u)) - params$b)
}

hess = function(u, params) {
    return(params$Q_beta)
}

# u_k for each partition
u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>% 
    group_by(leaf_id) %>% filter(l1_cost == min(l1_cost)) %>% 
    data.frame

# compute lambda_k for each of the u_k
# each lambda_k is stored column-wise, (D x K) matrix
lambda_k = apply(psi_df[,1:D], 1, grad, params = post)

# compute H_k for each of the u_k
H_k = Q_beta # same hessian for all partitions A_k

# compute m_k for each partition
m_k = mu_beta # same for all partitions A_k

## compute each of the terms in the exponential

log_terms = numeric(K)
G_k = numeric(K)
for (k in 1:K) {
    
    u_k = unname(unlist(psi_df[k,1:D]))
    diff_k = u_k - m_k
    
    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
    # G_k[k] = epmgp::pmvn(lb, ub, mu_beta, Q_beta_inv, log = TRUE)
    G_k[k] = log(TruncatedNormal::pmvnorm(m_k, Q_beta_inv, lb, ub)[1])
    
    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) - psi_df$psi_u[k] + 
        sum(lambda_k[,k] * u_k) - 
        0.5 * t(u_k) %*% H_k %*% u_k + 
        0.5 * t(m_k) %*% H_k %*% m_k + 
        G_k[k]
        
}

log_sum_exp(log_terms)


(LIL = lil(prior, post)) 










