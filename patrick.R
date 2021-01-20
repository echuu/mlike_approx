


# compute G_k for each partition
K = nrow(bounds)
G_k = numeric(K)
for (k in 1:K) {
    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist
    G_k[k] = epmgp::pmvn(lb, ub, mu_beta, Q_beta_inv, log = F)
    
}
G_k

library(TruncatedNormal)
p_0 = TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, lb, ub)
p_0[1]





