


# compute R(u)
# for a given u *in a specified partition*, compute psi(u) - psi_hat(u)
# where psi_hat(u) is the optimal psi value chosen in the previous layer

## start here
n_samps = 10

# for the partition learned from prev fitted tree, extract the partition id and
# the optimal value of psi for this partition
og_part = hml_approx$param_out %>% 
    dplyr::select(-c(psi_choice, logQ_cstar))

set.seed(1)
ss_part = fit_resid(og_part, D, n_samps, prior)
ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
fs_part = fit_resid(ts_part, D, n_samps / 2, prior)


# truth
(true_logml = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1]))

# original approx
hml_approx$const_vec
log_sum_exp(unlist(compute_expterms(ss_part, D)))
log_sum_exp(unlist(compute_expterms(ts_part, D)))
log_sum_exp(unlist(compute_expterms(fs_part, D)))


# repeat the same calculations above for G-many replications -------------------

G = 50
true_logml  = numeric(G)
orig_approx = numeric(G)
ss_approx   = numeric(G)
ts_approx   = numeric(G)

for (g in 1:G) {
    
    X   = matrix(rnorm(N * D), N, D)                # (N x D) design matrix
    eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))    # (N x 1) errors vector
    y   = X %*% beta + eps                          # (N x 1) response vector
    
    # compute posterior parameters ---------------------------------------------
    Q_beta     =  1 / sigmasq * (t(X) %*% X + tau * diag(1, D))
    Q_beta_inv =  solve(Q_beta)
    b          =  1 / sigmasq * t(X) %*% y
    mu_beta    =  Q_beta_inv %*% b
    
    # create prior, post objects to be passed into the hml algorithm 
    prior = list(y = y, X = X, sigmasq = sigmasq, tau = tau, N = N, D = D,
                 Q_beta = Q_beta, b = b)
    
    lil_0 = -0.5 * N * log(2 * pi) - 0.5 * (N + D) * log(sigmasq) + 
        0.5 * D * log(tau) - 0.5 * log_det(Q_beta) - 
        1 / (2 * sigmasq) * sum(y^2) + 0.5 * sum(b * mu_beta)
    
    true_logml[g] = lil_0 + D * log(2) + 
        log(TruncatedNormal::pmvnorm(mu_beta, Q_beta_inv, 
                                     lb = rep(0, D), ub = rep(Inf, D))[1])
    
    samples = data.frame(rtmvnorm(J, c(mu_beta), Q_beta_inv, 
                                  rep(0, D), rep(Inf, D)))
    
    u_df = preprocess(samples, D, prior)
    
    hml_approx = hml_const(1, D, u_df, J, prior)
    
    og_part = hml_approx$param_out %>% 
        dplyr::select(-c(psi_choice, logQ_cstar))
    
    ss_part = fit_resid(og_part, D, n_samps, prior)
    ts_part = fit_resid(ss_part, D, n_samps / 2, prior)
    
    orig_approx[g] = hml_approx$const_vec 
    
    
    ss_approx[g] = log_sum_exp(unlist(compute_expterms(ss_part, D)))
    ts_approx[g] = log_sum_exp(unlist(compute_expterms(ts_part, D)))
    
    # print(paste('iter ', g, '/', G, ' : ', 
    #             round(ss_approx[g], 4), ' (err: ', 
    #             round(abs(true_logml[g] - ss_approx[g]), 4), ', avg: ', 
    #             round(mean(ss_approx[1:g]), 4), '), ',
    #             round(ts_approx[g], 4), ' (err: ', 
    #             round(abs(true_logml[g] - ts_approx[g]), 4), ', avg: ', 
    #             round(mean(ts_approx[1:g]), 4), '), ', sep = ''))
    
    print(paste('iter ', g, '/', G, ' : ', 
                round(ss_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ss_approx[g]), 4), '), ',
                round(ts_approx[g], 4), ' (err: ', 
                round(abs(true_logml[g] - ts_approx[g]), 4), '), ', 
                sep = ''))
    
}

mean(true_logml)  # mean of the true log ML over all partitions
mean(orig_approx) # no re-partioning
mean(ss_approx)   # second stage 
mean(ts_approx)   # third stage

abs(mean(true_logml) - mean(orig_approx))
abs(mean(true_logml) - mean(ss_approx))
abs(mean(true_logml) - mean(ts_approx))







#### end of replications -------------------------------------------------------



part_0 = hml_obj$param_out %>% 
    dplyr::select(-c(psi_choice, psi_star, logQ_cstar))

part_set = part_0$leaf_id

orig_partition = hml_obj$param_out %>%
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs)) %>% 
    arrange(desc(perc))

K = length(part_set)

# initialize a list to store the vector containing the terms in exponential
# for each of the sub-partitions
# kth elmt is an s_k dim vector of terms that are to be exponentiated
# at the very end, all entries are unlisted and evaluated via log-sum-exp

exp_terms = vector("list", K) 
ck_star_list = vector("list", K)
ss_partition = vector("list", K) # store each of the hml_obj objects

perc_thresh = sort(orig_partition$perc, decreasing = T)

for (k in 1:K) {
    
    N_k_p = part_0$n_obs[k] * n_samps  # num of (re)samples to draw from part k
    part_k = part_0 %>%           # set of lower/upper bounds
        dplyr::filter(leaf_id == part_set[k]) %>%
        dplyr::select(-c(leaf_id, n_obs))
    
    # sample uniformly from each lower/upper bound pair to form a D-dim vector
    part_k_long = c(unlist(part_k)) %>% matrix(ncol = 2, byrow = T)
    
    resamp_k = Matrix_runif(N_k_p, lower = part_k_long[,1],
                            upper = part_k_long[,2]) %>% data.frame
    
    u_df_k = preprocess(resamp_k, D, prior) # N_k_p x (D_u + 1)
    
    c_k_approx = hml_const_mod(1, D, u_df_k, N_k_p, prior)
    
    ck_star_list[[k]] = c_k_approx$param_out %>%
        dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs)
    
    ss_partition[[k]] = c_k_approx
    
    exp_terms[[k]] = c_k_approx$const_approx
    
}





