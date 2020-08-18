
setwd("C:/Users/ericc/mlike_approx/algo")
source("setup.R")           

# load this LAST to overwrite def preprocess()
source("C:/Users/ericc/mlike_approx/mvn_ig/mvn_ig_helper.R") 




J = 1000          # number of MC samples per approximation
D = 20
N = 200 # for testing -- comment this line to perform ext. analysis
n_samps = 10

set.seed(123)
p       = D - 1            # dimension of beta
mu_beta = rep(0, p)        # prior mean for beta
V_beta  = diag(1, p)       # scaled precision matrix for betaV_beta
a_0     = 2 / 2            # shape param for sigmasq
b_0     = 1 / 2            # scale param 
beta    = sample(-10:10, p, replace = T)
sigmasq = 4                # true variance (1 x 1) 

I_p = diag(1, p)           # (p x p) identity matrix
I_N = diag(1, N)           # (N x N) identity matrix

X = matrix(rnorm(N * p), N, p) # (N x p) design matrix

eps = rnorm(N, mean = 0, sd = sqrt(sigmasq))

y = X %*% beta + eps # (N x 1) response vector
# ------------------------------------------------------------------


## compute posterior parameters ------------------------------------
V_beta_inv = solve(V_beta)
V_star_inv = t(X) %*% X + V_beta_inv

V_star  = solve(V_star_inv)                                # (p x p)
mu_star = V_star %*% (t(X) %*% y + V_beta_inv %*% mu_beta) # (p x 1)
a_n =  a_0 + N / 2 
b_n =  c(b_0 + 0.5 * (t(y) %*% y + 
                          t(mu_beta) %*% V_beta_inv %*% mu_beta - 
                          t(mu_star) %*% V_star_inv %*% mu_star))

# compute MLE estimates for mean, variance of the regression model
ybar = X %*% mu_star
sigmasq_mle = 1 / N * sum((y - ybar)^2)

# create prior, posterior objects
prior = list(V_beta = V_beta, 
             mu_beta = mu_beta, 
             a_0 = a_0, 
             b_0 = b_0,
             y = y, X = X,
             V_beta_inv = V_beta_inv)

# store posterior parameters
post  = list(V_star  =  V_star,
             mu_star =  mu_star,
             a_n     =  a_n,
             b_n     =  b_n,
             V_star_inv = V_star_inv)


# sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
# beta_post = matrix(0, J, p)
# for (j in 1:J) {
#     beta_post[j,] = rmvnorm(1, mean = mu_star, sigma = sigmasq_post[j] * V_star)
# }
sample_beta = function(s2, post) {
    rmvnorm(1, mean = post$mu_star, sigma = s2 * post$V_star)
}

sigmasq_post = MCMCpack::rinvgamma(J, shape = a_n, scale = b_n)
beta_mat = t(sapply(sigmasq_post, sample_beta, post = post))
u_samp = data.frame(beta_mat, sigmasq_post)
u_df = preprocess(u_samp, D, prior)

# refactored version
hml_approx = hml_const(1, D, u_df, J, prior) 
hml_approx$param_out %>% 
    dplyr::select(leaf_id, psi_choice, psi_star, logQ_cstar, n_obs) %>% 
    dplyr::mutate(perc = n_obs / sum(n_obs))

hml_approx$const_vec      # -256.761
lil(y, X, prior, post)    # -256.7659

og_part = hml_approx$param_out %>%
    dplyr::select(-c(psi_choice, logQ_cstar))
ss_part = fit_resid(og_part, D, n_samps, prior)
ts_part = fit_resid(ss_part, D, n_samps / 2, prior)

x1 = hml_approx$const_vec
x2 = log_sum_exp(unlist(compute_expterms(ss_part, D)))
x3 = log_sum_exp(unlist(compute_expterms(ts_part, D)))

(approx = c(x1, x2, x3))


### check containment for each point
matrix_part = function(row_part) {
    row_part %>% matrix(ncol = 2, byrow = TRUE)
}
check_member = function(u, A_k) {
    all(A_k[,1] <= u) && all(A_k[,2] >= u)
}
## given a point, determine which of the k partitions it belongs in
query_partition = function(u, part_list) {
    ind = (1:length(part_list))[(sapply(part_list, check_member, u = u))]
    if (length(ind) == 0) {
        # print("looking in old partition")
        # ind = (1:length(og_part_list))[(sapply(og_part_list, check_member, u = u))]
        ind = NA
    }
    return(ind)
}


############# start weighting here ---------------------------------------------


# extract posterior samples without Psi(u) evaluation
u_sub = u_df %>% dplyr::select(-psi_u)

############ process stage 1 ---------------------------------------------------

# extract only partition-defining columns
part1 = og_part %>% dplyr::select(-c('leaf_id', 'psi_star', 'n_obs'))

# extract the optimal value for each partition
psi_star_1 = og_part$psi_star

# convert partitions stored as K rows into a list of K partitions
part1_list = lapply(split(part1,seq(NROW(part1))), matrix_part)
part1_id = apply(u_sub, 1, query_partition, part_list = part1_list)
part1_id %>% head

u_df_wt = u_df %>% dplyr::mutate(psi_1 = psi_star_1[part1_id])

############ process stage 2 ---------------------------------------------------

# extract only partition-defining columns
part2 = ss_part %>% dplyr::select(-c('leaf_id', 'psi_star', 'n_obs'))

# extract the optimal value for each partition
psi_star_2 = ss_part$psi_star

# convert partitions stored as K rows into a list of K partitions
part2_list = lapply(split(part2, seq(NROW(part2))), matrix_part)

part2_id = apply(u_sub, 1, query_partition, part_list = part2_list)

u_df_wt = u_df_wt %>% dplyr::mutate(psi_2 = psi_star_2[part2_id])

# replace missing values with psi_1 values (this will be done for stage 3 too)
u_df_wt = u_df_wt %>% dplyr::mutate(psi_2 = ifelse(is.na(psi_2), psi_1, psi_2))

u_df_wt$psi_2 %>% is.na %>% sum

############ process stage 3 ---------------------------------------------------

# extract only partition-defining columns
part3 = ts_part %>% dplyr::select(-c('leaf_id', 'psi_star', 'n_obs'))

# extract the optimal value for each partition
psi_star_3 = ts_part$psi_star

# convert partitions stored as K rows into a list of K partitions
part3_list = lapply(split(part3, seq(NROW(part3))), matrix_part)

part3_id = apply(u_sub, 1, query_partition, part_list = part3_list)

sum(is.na(part3_id))

u_df_wt = u_df_wt %>% dplyr::mutate(psi_3 = psi_star_3[part3_id])

# replace missing values with psi_1 values (this will be done for stage 3 too)
u_df_wt = u_df_wt %>% dplyr::mutate(psi_3 = ifelse(is.na(psi_3), psi_2, psi_3))

u_df_wt$psi_3 %>% is.na %>% sum

u_df_wt %>% dplyr::select(psi_u, psi_1, psi_2, psi_3) %>% head()

############ compute weights for each stage ------------------------------------
library(MLmetrics)
# (d1 = MSE(u_df_wt$psi_u, u_df_wt$psi_1))
# (d2 = MSE(u_df_wt$psi_u, u_df_wt$psi_2))
# (d3 = MSE(u_df_wt$psi_u, u_df_wt$psi_3))

(d1 = MAE(u_df_wt$psi_u, u_df_wt$psi_1))
(d2 = MAE(u_df_wt$psi_u, u_df_wt$psi_2))
(d3 = MAE(u_df_wt$psi_u, u_df_wt$psi_3))

approx # 1/2/3 stage approx
wt1 = exp(-d1 / D)
wt2 = exp(-d2 / D)
wt3 = exp(-d3 / D)
wt_norm = wt1 + wt2 + wt3

(wt_vec = c(wt1, wt2, wt3) / wt_norm)

sum(wt_vec * approx)    # exp weights
mean(approx)            # average
lil(y, X, prior, post)  # truth

# matrix_part(part1[1,])
# apply(part1, 1, matrix_part)
# j = 1
# 
# u1 = u_sub[j,]
# A1 = part1[k,] %>% matrix(ncol = 2, byrow = TRUE)
# 
# k = 2
# K = length(part1_list)
# membership = as.logical(numeric(K))
# for (k in 1:K) {
#     membership[k] = all(part1_list[[k]][,1] <= u1) && 
#         all(part1_list[[k]][,2] >= u1)
# }

















