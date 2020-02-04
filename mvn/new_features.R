

DELL_PATH = "C:/Users/chuu/mlike_approx"
setwd(DELL_PATH)

source("partition/partition.R")  # load partition extraction functions
source("extractPartition.R")     # extractPartition() function
source("hybrid_approx.R")        # load main algorithm functions
source("mvn/mvn_helper.R")       # load psi(), lambda()


x11()

library(pracma)

fun = function(x, y) exp(-1/2 * (x^2 + y^2))



# model settings ---------------------------------------------------------------

D  = 2
Sigma = diag(1, D)
Sigma_inv = solve(Sigma)
prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv)
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500

# ------------------------------------------------------------------------------


# (1) demonstrate use of the algorithm for a given D, N ------------------------

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

# ------------------------------------------------------------------------------


# (2) extract details of the approximation by partition ------------------------


psi_fun = fun # move this back into the function call 
mvn_diag = approx_lil_diag(D, u_df_full, prior) # hybrid_approx.R

mvn_diag$logZ_numer
mvn_diag$logZ_taylor1
mvn_diag$lozZ_taylor2
mvn_diag$partition_info
mvn_diag$taylor2_integral
mvn_diag$verbose_partition

partition_info = mvn_diag$partition_info %>% 
    mutate(numer = round(numer, 4), taylor1 = round(taylor1, 4), 
           lambda1 = round(lambda1, 5), lambda2 = round(lambda2, 5), 
           taylor2 = round(taylor2, 4), e_ck_2 = round(e_ck_2, 4))

partition_info

# write.csv(partition_info, "partition_info_mvn.csv", 
#          row.names = F)

plotPartition(u_df_full, mvn_diag$param_out)


# ------------------------------------------------------------------------------

# TODO: residual analysis for each partition, correct the median residual 
#       to the MIN residual

# (1) for each partition, compute the (sum of) residuals using both of the
#     approximation schemes - constant approx, 1st order taylor. store these
#     as columns in the 'partition_info' dataframe. next, these will be used
#     to determine which of the two approximations to use (AB suspects that
#     because of the formulation of constant approximation, the constant approx
#     will always beat out the 1st order taylor)


#### the code below will go into a NEW function that separate from the main
#### algorithm -- for now, it will be used as a post-process approximation
#### method to improve the algorithm 
#### eventually, this should be implemented within the algorithm, but for
#### illustrative purposes, we still want BOTH constant and taylor
#### approximations so we can see how this new way of approximations performs

# extract the rpart object that's used in the main algorithm
mvn_rpart = mvn_diag$u_rpart
param_out = mvn_diag$param_out

star_ind = grep("_star", names(param_out)) # columns of u_k_star for each k

# appends the leaf id for each observation as a column
u_df = u_df_full %>% mutate(leaf_id = mvn_rpart$where, 
                            const_approx = 0,  const_resid = 0, 
                            taylor_approx = 0, taylor_resid = 0)

partition_id = mvn_rpart$where %>% unique
n_partitions = length(table(u_df$leaf_id))
for (i in 1:n_partitions) {
    
    k = partition_id[i]
    
    u_k_star = param_out %>% filter(leaf_id == k) %>% select(star_ind) %>% 
        unname %>% unlist
    
    # compute constant approximation for psi
    u_df[u_df$leaf_id == k,]$const_approx = psi(u_k_star, prior) %>% c()
    
    # compute order 1 taylor approximation for psi
    # 
    diff_k = sweep(u_df %>% filter(leaf_id == k) %>% select(u1, u2), 2, 
                   FUN = '+', -u_k_star)             
                                        
    u_df[u_df$leaf_id == k,]$taylor_approx = c(psi(u_k_star, prior)) + 
        as.matrix(diff_k) %*% lambda(u_k_star, prior)
    
    u_df = u_df %>% mutate(const_resid  = psi_u - const_approx,
                           taylor_resid = psi_u - taylor_approx)
    
}

error_df = data.frame(leaf_id = partition_id, 
                      const_sq_error = 0, taylor_sq_error = 0)

for (i in 1:n_partitions) {
    
    k = partition_id[i]
    
    # const_sq_error
    sse_const  = sum(u_df[u_df$leaf_id == k,]$const_resid^2)
    sse_taylor = sum(u_df[u_df$leaf_id == k,]$taylor_resid^2)
    
    # for each partition, calulcate sum of residuals for const and taylor approx
    error_df = error_df %>% 
        mutate(const_sq_error = replace(const_sq_error, partition_id == k, 
                                        sse_const),
               taylor_sq_error = replace(taylor_sq_error, partition_id == k,
                                         sse_taylor))
}


partition_approx = mvn_diag$verbose_partition %>% 
    select(leaf_id, numer, taylor1, taylor2)

partition_approx = merge(partition_approx, error_df, by = "leaf_id")

# extract leaf id for which we use the taylor approximation
taylor_index = error_df %>% filter(taylor_sq_error < const_sq_error) %>% 
    select(leaf_id) %>% unlist %>% unname

# extract leaf id for which we use the constant approximation
const_index = error_df %>% filter(taylor_sq_error >= const_sq_error) %>% 
    select(leaf_id) %>% unlist %>% unname


mvn_diag$partition_info %>% select(numer) %>% sum %>% log

const_contribution = mvn_diag$verbose_partition %>% 
    filter(leaf_id %in% const_index) %>% 
    select(taylor1)

taylor_contribution = mvn_diag$verbose_partition %>% 
    filter(leaf_id %in% taylor_index) %>% 
    select(taylor2)

hybrid_approx = log(sum(const_contribution, taylor_contribution)) # 2.113237

#### end of additional algorithm code ------------------------------------------

# compare with constant approximation and taylor approximation

mvn_diag$logZ_taylor1 # 2.102216

mvn_diag$lozZ_taylor2 # 3.311019


# (2) in current implementation, we use the point whose residual is MEDIAN, but
#     this doesn't make much intuitive sense - instead, use point whose
#     residual is MINIMUM

# (3) build model-based tree with the following specification 
#     y ~ u1, u2, ... , uD | u1, u2, ... , uD

# (4) run the algorithm for the MCMC model - recall, D = 2 since this is 
#     dimension of u, but there are 20 parameters in the mixture density



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





