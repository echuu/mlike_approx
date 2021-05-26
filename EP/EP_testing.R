
D = D_u
## (2) fit the regression tree via rpart()
u_rpart = rpart::rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)

#### extension starts here -------------------------------------------------

### (1) find global mean
MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
u_star = u_0

u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_star)
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

# k = 2
k = 1
for (k in 1:K) {

    u_k = unname(unlist(psi_df[k,1:D]))

    H_k = pracma::hessian(psi, u_k, params = params)
    # H_k = hess(u_k, params)
    H_k_inv = chol2inv(chol(H_k))

    lambda_k = pracma::grad(psi, u_k, params = params)
    # lambda_k = grad(u_k, params = params)
    b_k = H_k %*% u_k - lambda_k
    m_k = H_k_inv %*% b_k

    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

    epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])
    source("R/epmgp.R")
    hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub)$logZ


    # do EP w/ modified step for updating site parameters and compare with the
    # estimates given above from standard EP, and botev tilted


    # G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # G_k[k] = ep_step(lb, ub, m_k, H_k_inv)
    G_k[k] = epmgp_stable(m_k, H_k_inv, b_k, lb, ub)$logZ
    # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
        psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k +
        0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]

}

G_k
log_sum_exp(log_terms)



(LIL = logmarginal(Y, testG, b, V, S))
