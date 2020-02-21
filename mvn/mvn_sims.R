


library("mvtnorm")           # for draws from multivariate normal

# path for lenovo
LEN_PATH  = "C:/Users/ericc/mlike_approx"
setwd(LEN_PATH)

# path for dell
# DELL_PATH = "C:/Users/chuu/mlike_approx"
# setwd(DELL_PATH)

# files must be loaded in this order, since *_helper.R file will sometimes
# overwrite a functino defined in hybrid_approx.R depending on the example

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("mvn/mvn_helper.R")              # load psi(), lambda() function
source('extractPartition.R')            # load extractPartition() function

D  = 10
Omega = diag(1, D)
N  = 200
Sigma = D / N * Omega 
Sigma_inv = solve(Sigma)
mu_0 = rep(0, D)
prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, mu_0 = mu_0)

# closed form of the normalizing constant
(logZ_mvn = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma))


B        = 100          # number of batch estimates
J        = 100          # number of MC samples
N_approx = 1            # number of estimates to compute per iteration b

set.seed(1)
u_samps = rmvnorm(J * B, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df = preprocess(u_samps, D, prior)
u_df %>% head

hml_approx = hml(N_approx, D, u_df, J, prior)
hml_approx$hybrid_vec


# ------------------------------------------------------------------------------


lil_hme  = numeric(B) # store the harmonic mean approximate for each batch
prop_const = numeric(B) # proportion of partitions that use constant approx 


for (b_i in 1:B) {
    
    if (b_i %% 20 == 0) print(paste("iter: ", b_i, "/", B, sep = ""))
    
    start = ((b_i - 1) * J + 1)
    end = start + J - 1
    
    u_df_b = u_df[start:end,]
    
    # (1) compute hybrid app
    hml_approx = hml(N_approx, D, u_df_b, J, prior) 
    lil_hml[b_i] = hml_approx$hybrid_vec
    
    prop_const[b_i] = hml_approx$n_const / nrow(hml_approx$verbose_partition)
    
    # (2) harmonic mean estimator
    # lil_hme[b_i] = hme(u_df, prior, J, D, N)
    
} # end of simulation outer loop


## plot results
# hme_df = data.frame(mcmc = 1:B, hml = lil_hml, hme = lil_hme, came = lil_came)
hme_df = data.frame(mcmc = 1:B, hml = lil_hml)

# mean average error (AE, true - estimated)
(MAE = round(mean(logZ_mvn - lil_hml), 3))

# root mean squared error (RMSE)
(RMSE = round(sqrt(mean((logZ_mvn - lil_hml)^2)), 3))

ggplot(hme_df, aes(x = mcmc, y = lil_hml)) + geom_point(col = 'blue') +
    geom_hline(aes(yintercept = logZ_mvn), linetype = 'dashed', size = 1.3) +
    labs(x = 'iter', y = 'logML',
         title = paste("logZ = ", round(logZ_mvn, 3), 
                       ", D = ", D, ", J = ", J, 
                       " (MAE = ", MAE, ", RMSE = ", RMSE, ")", sep = ''))

prop_const # prop of constant approximations used for each partition
mean(prop_const)
























