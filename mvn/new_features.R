
source("partition/partition.R")  # load partition extraction functions
source("extractPartition.R")
source("hybrid_approx.R")        # load main algorithm functions
source("mvn/mvn_helper.R")       # load psi(), lambda()


x11()

library(pracma)

fun = function(x, y) exp(-1/2 * (x^2 + y^2))

# taylor_int = function(u1, u2) { # 2nd term of the taylor expansion
#     exp(-(lambda1 * (u1 - u_k_star[1]) + lambda2 * (u2 - u_k_star[2])))
# }

# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

D  = 2
Sigma = diag(1, D)
Sigma_inv = solve(Sigma)


prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv)
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500

set.seed(1)
J = 5000
N_approx = 1
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew) # code from mvn.R should match this

# compare approximation with the true integral
result = pracma::integral2(fun, -100, 100, -100, 100, reltol = 1e-50)
log(result$Q) # 1.837877

psi_fun = fun

mvn_diag = approx_lil_diag(D, u_df_full, prior)

mvn_diag$logZ_numer
mvn_diag$logZ_taylor1
mvn_diag$lozZ_taylor2
mvn_diag$partition_info
mvn_diag$param_out
mvn_diag$taylor2_integral
mvn_diag$verbose_partition

partition_info = mvn_diag$partition_info %>% 
    mutate(numer = round(numer, 4), taylor1 = round(taylor1, 4), 
           lambda1 = round(lambda1, 5), lambda2 = round(lambda2, 5), 
           taylor2 = round(taylor2, 4), e_ck_2 = round(e_ck_2, 4))


partition_info

write.csv(partition_info, "partition_info_mvn.csv", 
          row.names = F)



plotPartition(u_df_full, mvn_diag$param_out)



# ------------------------------------------------------------------------------

library(tree)
u_tree = tree(psi_u ~ ., u_df_full)
plot(u_tree)
text(u_tree, cex = 0.8)

plot(u_df_full[,1], u_df_full[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'u1', ylab = 'u2', main = '')
partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


#### write out the new features in the chunk below

u_rpart = rpart(psi_u ~ ., u_df_full)
#plot(u_rpart)
#text(u_rpart, cex = 0.7)

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df_full[,d])
    param_d_max = max(u_df_full[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)  # partition.R
# organize all data into single data frame --> ready for approximation
param_out = u_star(u_rpart, u_df_full, u_partition, D)
K = nrow(u_partition)     # number of partitions 


# 2nd term of the taylor expansion
taylor_int = function(u1, u2) {
    exp(-(lambda1 * (u1 - u_k_star[1]) + lambda2 * (u2 - u_k_star[2])))
}



# newly declared storage
partition_integral = numeric(K) # store the numerical integral of each partition
area_k             = numeric(K) # store the area of each partition

taylor1_approx = numeric(K)     # store approx that uses 1-term taylor
taylor2_approx = numeric(K)     # store approx taht uses 2-term taylor

e_ck_1   = numeric(K)           # store first constant in taylor approx
e_ck_2   = numeric(K)           # store second constant in taylor approx


lambda_k = data.frame(matrix(NA, K, D)) # store gradient evaluated at u_k_star
taylor2_integral = numeric(K)       # store integral of e^(-lambda_k'u) over A_k


star_ind = grep("_star", names(param_out)) # columns of u_k_star for each k

for (k in 1:K) {
    
    u1_b = param_out[k, 6]  # upper bound of u1
    u1_a = param_out[k, 5]  # lower bound of u1
    
    u2_b = param_out[k, 8]  # upper bound of u2
    u2_a = param_out[k, 7]  # lower bound of u2
    
    # (0) true value of the integral over the k-th partition
    result = integral2(fun, u1_a, u1_b, u2_a, u2_b, reltol = 1e-50)
    partition_integral[k] = result$Q
    
    # --------------------------------------------------------------------------
    
    # for each parttion, A_k, keep track of the following quantities: 
    #     (1) one term taylor approximation
    #         (1.1) e_ck_1 = e^(-psi(u_kstar))
    #         (1.2) simple approx: e_ck_0 * (Area of A_k partition)
    #     (2) two term taylor approximation
    #         (2.1) lambda_k = gradient evaluated at u_k_star
    #         (2.2) e_ck_2 = e^(lambda_k'u_k_star) = 2nd (constant) term
    #         (2.3) taylor2_integral = integral of e^(-lambda_k'u) over A_k
    #         (2.4) taylor2_approx = e_ck_1 * e_ck_2 * taylor2_integral
    
    # (1) compute integral via one term taylor approximation
    u_k_star = param_out[k, star_ind] %>% unlist %>% unname
    
    # e_ck[k] = exp(-psi(u_k_star, prior))
    area_k[k] = (u1_b - u1_a) * (u2_b - u2_a)
    
    e_ck_1[k] = exp(-psi(u_k_star, prior))
    
    taylor1_approx[k] = e_ck_1[k] * area_k[k]
    
    # (2) compute integral via two term taylor approximation
    
    lambda_k[k,] = lambda(u_k_star, prior)
    lambda1  = lambda_k[k, 1]
    lambda2  = lambda_k[k, 2]
    
    
    # (2.2) e_ck_2 = e^(lambda_k'u_k_star) = 2nd (constant) term
    e_ck_2[k] = exp(sum(lambda_k[k,] * u_k_star))
    
    # (2.3) taylor2_integral = integral of e^(-lambda_k'u) over A_k
    # normally this is D-iteration loop that calculates each 1-d integral
    taylor2_integral[k]  = - 1 / lambda1 * - 1 / lambda2 * 
        (exp(-lambda1 * u1_b) - exp(-lambda1 * u1_a)) * 
        (exp(-lambda2 * u2_b) - exp(-lambda2 * u2_a))
    
    taylor2_approx[k] = e_ck_1[k] * e_ck_2[k] * taylor2_integral[k]
    
    # (3) compute nmerical integral over these regions
    # order1_numer[k]  = integral2(taylor_int, u1_a, u1_b, u2_a, u2_b, 
    #                              reltol = 1e-50)$Q
    # taylor2_numer[k] = e_ck_1[k] * order1_numer[k]
    
}


log(sum(partition_integral))  # numerically calculated integrals over partition
# log(sum(taylor2_numer))     # both match for case when D = 2, N not used
log(sum(taylor2_approx))      # approximation that comes out of approx_lil()


log(sum(taylor1_approx))
log(sum(taylor1_approx[7:14]) + sum(taylor2_approx[1:6]))
log(sum(taylor1_approx[4:14]) + sum(taylor2_approx[1:4]))


names(lambda_k) = c("lambda1", "lambda2")

all_integrals = cbind(numer = partition_integral, taylor1 = taylor1_approx) %>% 
    cbind(taylor2 = taylor2_approx)

diagnostics = all_integrals %>% cbind(lambda_k) %>% cbind(e_ck_2)

diagnostics

# look at the integrals (and their approximations) over each of the partitions.
# note, this is sorted by the percent membership in each of the leaf nodes.
# in general, larger membership -> smaller psi_hat -> larger integral

(param_out %>% mutate(perc_mem = n_obs / sum(n_obs)))[,c(1:4,10)] %>% 
    select(perc_mem, psi_hat, u1_star, u2_star) %>% cbind(all_integrals) %>% 
    arrange(desc(perc_mem))

(param_out %>% mutate(perc_mem = n_obs / sum(n_obs))) %>% 
    select(perc_mem, psi_hat, u1_star, u2_star) %>% cbind(diagnostics) %>% 
    arrange(desc(perc_mem)) %>% round(4)





#### column description:
# perc_mem : percent membership
# psi_hat  : predicted value from the tree
# u1_star  : "representative point" for 1st component
# u2_star  : "representative point" for 2nd component
# numer    : numerically computed integral over each partition
# taylor1  : integral of the 1-term taylor approximation over each partition
# taylor2  : integral of the 2-term taylor approximation over each partition
# lambda1  : 1st component of gradient evaluated at u_k_star
# lambda2  : 2nd component of gradient evaluated at u_k_star
# e_ck_2   : 2nd constant in the 2-term taylor approximation





