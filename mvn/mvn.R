



library("mvtnorm")           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation

# LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
# setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("mvn/mvn_helper.R")  # load psi(), lambda()



D  = 2
N  = 500
Sigma = D / N * diag(1, D)
Sigma_inv = solve(Sigma)
mu_0 = rep(0, D)

prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, mu_0 = mu_0)

D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500


set.seed(1)
J = 5000
N_approx = 10
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew)


u_tree = tree(psi_u ~ ., u_df_full)



# step through of algorithm to see values of zhat ------------------------------

plot(u_tree)
text(u_tree, cex = 0.8)

plot(u_df_full[,1], u_df_full[,2], pch = 20, cex = 0.8, col = "cyan",
     xlab = 'u1', ylab = 'u2', main = '')
partition.tree(u_tree, add = TRUE, cex = 0.8, ordvars = c("u1", "u2"))


# start code here for mvn_troubleshoot.R ---------------------------------------

u_df = u_df_full

u_rpart = rpart(psi_u ~ ., u_df)
#plot(u_rpart)
#text(u_rpart, cex = 0.7)



# (3.1) obtain the (data-defined) support for each of the parameters
param_support = matrix(NA, D, 2) # store the parameter supports row-wise

for (d in 1:D) {
    param_d_min = min(u_df[,d])
    param_d_max = max(u_df[,d])
    
    param_support[d,] = c(param_d_min, param_d_max)
}

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)  # partition.R

# organize all data into single data frame --> ready for approximation
param_out = u_star(u_rpart, u_df, u_partition, D)

n_partitions = nrow(u_partition)     # numebr of partitions 
c_k          = numeric(n_partitions) # constant term for k-th partition
zhat         = numeric(n_partitions) # integral over k-th partition

# (4) compute closed form integral over each partition
lambda_mat = matrix(NA, n_partitions, D)

for (k in 1:n_partitions) {
    
    # extract "representative point" of the k-th partition
    star_ind = grep("_star", names(param_out))
    u = param_out[k, star_ind] %>% unlist %>% unname
    
    l_k = lambda(u, prior)       # (D x 1) 
    
    # evaluate e^c_k = e^{psi(u_star)}
    c_k[k] = exp(-psi(u, prior)) # (1 x 1)
    c_k[k] = exp(-psi(u, prior) + sum(l_k * u)) 
    
    # compute lambda_k : gradient of psi, evaluated at u_star
    lambda_mat[k,] = l_k
    
    # store each component of the D-dim integral 
    integral_d = numeric(D)      # (D x 1)
    
    for (d in 1:D) {
        
        # find column id of the first lower bound
        col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
        col_id_ub = col_id_lb + 1
        
        upper_kd = param_out[k, col_id_ub]
        lower_kd = param_out[k, col_id_lb]
        
        # d-th integral computed in closed form
        integral_d[d] = - 1 / l_k[d] * 
             (exp(-l_k[d] * upper_kd) - exp(-l_k[d] * lower_kd)) 
        
        # area of the partition
        # integral_d[d] = (param_out[k, col_id_ub] - param_out[k, col_id_lb])
        
    } # end of loop computing each of 1-dim integrals
    
    # print(l_k)
    
    # compute the D-dim integral (product of D 1-dim integrals)
    zhat[k] = prod(c_k[k], integral_d)
    # zhat[k] = prod(integral_d)
    
} # end of for loop over the K partitions

log(sum(zhat))






cbind(param_out[,1:4], zhat) %>% cbind(lambda_mat)
integral_d
zhat










# move stuff below (w/ necessary params defined above) into asymptotics.R file
# ------------------------------------------------------------------------------
N_vec_log = seq(6, 12, 0.05)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
logZ_0    = numeric(length(N_vec))
for (i in 1:length(N_vec)) {
    
    
    N = N_vec[i]
    Sigma = D / N * diag(1, D)
    Sigma_inv = solve(Sigma)
    
    logZ_0[i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma)
    
    #print(paste("iter = ", i, 
    #            " -- Calculating LIL for D = ", D, ", N = ", N, sep = ''))
}


lil_df = data.frame(logZ_0 = logZ_0, logn = log(N_vec))

formula1 = y ~ x
ggplot(lil_df, aes(logn, logZ_0)) + geom_point() + 
    labs(title = "") + 
    geom_smooth(method = lm, se = T, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ_0~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16)





