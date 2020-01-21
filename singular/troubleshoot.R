

install.packages("pracma")
library(pracma)
fun <- function(x, y) exp(-n*x^2*y^4)
logn = seq(1,20,by=0.1)
N = exp(logn)
lML = rep(0,length(N))

i = 1;
n = 100
n = 1000
for (n in N){
    result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
    lML[i] = log(result$Q)
    i = i + 1;
}
plot(logn,lML)
lm(lML ~ logn)


n = 100
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -0.724

n = 1000
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014


N = 1000
D = 2                               # dimension of parameter

gamma_dat = list(N = N)

# should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  'gamma_sample.stan', 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   seed    =  stan_seed,
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)          

u_df_N = preprocess(gamma_fit_N, D, N)

rpart_sing = rpart(psi_u ~ ., u_df_N, cp = 0.001)
plot(rpart_sing)
text(rpart_sing, size = 0.5)

approx_test = approx_lil_stan(N_approx, D, N, u_df_N, J)

approx_test

mean(approx_test)

# ------------------------------------------------------------------------------

u_df_N %>% head

u_rpart = rpart(psi_u ~ ., u_df_in)

rpart_obj = u_rpart


u_rpart = sing_rpart

 
plot(u_df_N[,1], u_df_N[,2])

plot(sing_rpart)
text(sing_rpart, size = 0.5)



 # -----------------------------------------------------------------------------

# STAN SETTINGS ----------------------------------------------------------------
J         = 1000         # number of MC samples per approximation
N_approx  = 10           # number of approximations
burn_in   = 2000         # number of burn in draws
n_chains  = 4            # number of markov chains to run
stan_seed = 123          # seed

J_iter = 1 / n_chains * N_approx * J + burn_in 
# ------------------------------------------------------------------------------

# GLOBAL MODEL SETTINGS --------------------------------------------------------
D = 2                    # dimension of parameter

# one run of the algorithm -----------------------------------------------------
set.seed(1)
N = 1000

gamma_dat = list(N = N) # for STAN sampler
prior     = list(N = N) # for evaluation of psi, lambda

# (1) generate posterior samples -- should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  'singular/gamma_sample.stan', 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   seed    =  stan_seed,
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)           

u_samp = rstan::extract(gamma_fit_N, pars = c("u"), permuted = TRUE)
u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2

# (2) evaluate posterior samples using psi(u)
u_df_N = preprocess(u_post, D, prior)

u_df = u_df_N

u_rpart = rpart(psi_u ~ ., u_df)
plot(u_rpart)
text(u_rpart, cex = 0.7)


# plot the partition
library(tree)

u_tree = tree(psi_u ~ ., u_df_N)

plot(u_tree)
text(u_tree, cex = 0.8)
partition.tree(u_tree, cex = 1, ordvars = c("u1", "u2"))

plot(u_df_N[,1], u_df_N[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'u1', ylab = 'u2', main = '')
partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


# ------------------------------------------------------------------------------

# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition
u_partition = paramPartition(u_rpart, param_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(u_rpart, u_df, u_partition, D)

n_partitions = nrow(u_partition)     # numebr of partitions 
c_k          = numeric(n_partitions) # constant term for k-th partition
zhat         = numeric(n_partitions) # integral over k-th partition

# (4) compute closed form integral over each partition
for (k in 1:n_partitions) {
    
    # extract "representative point" of the k-th partition
    star_ind = grep("_star", names(param_out))
    u = param_out[k, star_ind] %>% unlist %>% unname
    
    # evaluate e^c_k = e^{psi(u_star)}
    c_k[k] = exp(-psi(u, prior)) # (1 x 1)
    
    # compute lambda_k : gradient of psi, evaluated at u_star
    l_k = lambda(u, prior)       # (D x 1) 
    
    # store each component of the D-dim integral 
    integral_d = numeric(D)      # (D x 1)
    
    for (d in 1:D) {
        
        # find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
        col_id_ub = col_id_lb + 1
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
            exp(- l_k[d] * (param_out[k, col_id_ub] - 
                                param_out[k, col_id_lb]))        
    } # end of loop computing each of 1-dim integrals
    
    print(integral_d)
    
    # compute the D-dim integral (product of D 1-dim integrals)
    zhat[k] = prod(c_k[k], integral_d)
    
} # end of for loop over the K partitions

# store the log integral \approx log marginal likelihood
log(sum(zhat))

# ------------------------------------------------------------------------------

n = 1000
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014


cbind(param_out[,1:4], zhat)


# 2-d contour plot of the parameter space
ggplot(u_df_N, aes(u1, u2)) + geom_density_2d()

ggplot(u_df_N, aes(u1, u2)) + geom_point()























