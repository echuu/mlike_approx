

# simulate data ------------------------------------------------------

N_vec = c(50, 100, 200, 500, 750, 1000, 2000, 5000, 8000, 10000)
J = 2000
D = 2
num_sims = length(N_vec)

u_mat = matrix(NA, J, num_sims * 2) # (J x 2 * num_sums)

set.seed(1)
# simulate all data first --> save 
for (i in 1:num_sims) {
    
    N = N_vec[i]
    
    print(paste("Generating data for N =", N))
    
    gamma_dat = list(N = N)
    
    # uncomment to refit model
    gamma_fit = stan(file = 'gamma_sample.stan', data = gamma_dat)
    
    u_samp = extract(gamma_fit, pars = c("u"), permuted = TRUE)
    
    u_post = u_samp$u[1:J,] %>% data.frame() # (J x 2) : post sample stored row-wise
    
    col_id = 2 * (i - 1) + 1
    
    u_mat[, col_id]     = u_post[, 1]
    u_mat[, col_id + 1] = u_post[, 2]
    
}


write.csv(u_mat, "gamma_samples.csv", row.names = F)
test_read = read.csv("gamma_samples.csv")
dim(u_mat)


gamma_dat = list(N = N)

# uncomment to refit model
gamma_fit = stan(file = 'gamma_sample.stan', data = gamma_dat,
                 seed = 1,)

u_samp = extract(gamma_fit, pars = c("u"), permuted = TRUE)

u_post = u_samp$u[1:J,] %>% data.frame() # (J x 2) : post sample stored row-wise

psi_u = apply(u_post, 1, psi, N = N)     # (J x 1) : negative log posterior

u_df_names = character(D + 1)
for (d in 1:D) {
    u_df_names[d] = paste("u", d, sep = '')
}
u_df_names[D + 1] = "psi_u"

# populate u_df
u_df = cbind(u_post, psi_u) # J x (D + 1)

# rename columns (needed since these are referenced explicitly in partition.R)
names(u_df) = u_df_names

u_rpart = rpart(psi_u ~ ., u_df)

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition --- moment of truth!!
u_partition = paramPartition(u_rpart, param_support)  # partition.R

head(u_partition)

# (3.3) organize all data into single data frame (see partition.R for format)

# extracts u_star, representative point of each partition (u_star \in R^D)
# psi_hat, leaf_id, u1_star, u2_star, ... , uD_star, 
#                   u1_lb, u1_ub, ...uD_lb, uD_ub
param_out = u_star(u_rpart, u_df, u_partition, D)

head(param_out)

# test lambda 
# k = 2
# u_test = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
# exp(-psi(u_test, N))     # (1 x 1) -- c_k[k] calculation
# lambda(u_test, N) # (D x 1) -- lambda(u_star) calculation


## (4) begin main algorithm 
n_partitions = nrow(u_partition)
c_k = numeric(n_partitions)
zhat = numeric(n_partitions)

for (k in 1:n_partitions) {
    
    # u_star_k = (mu_k, sigma_sq_k)
    # c_k[k] = exp(-psi(param_out[k,]$u1_star, 
    #                   param_out[k,]$u2_star,
    #                   y, m_0, w_0, r_0, s_0)) # (1 x 1)
    
    # TODO: avoid explicit definition of each representative point
    u = c(param_out[k,]$u1_star, param_out[k,]$u2_star)
    
    c_k[k] = exp(-psi(u, N)) # (1 x 1)
    
    l_k = lambda(u, N)
    
    integral_d = numeric(D) # store each component of the D-dim integral 
    
    # nothing to refactor in this loop (i think?) since we're just iterating
    # thru each of the integrals and computing an exponential term
    for (d in 1:D) {
        
        # verify these -- these need to be recalculated if the form of param_out
        # changes (if columns get shuffled)
        
        # col id will change for D > 2
        # DONE: generalize this better so there's less obscure calculation
        # col_id_lb = 5 + 2 * (d - 1)
        # col_id_ub = col_id_lb + 1
        
        # updated 1/14: find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out))
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
            exp(- l_k[d] * (param_out[k, col_id_ub] - param_out[k, col_id_lb]))        
        
    }
    
    zhat[k] = prod(c_k[k], integral_d)
}

log(sum(zhat))







