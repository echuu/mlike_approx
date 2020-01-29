
## mvn.R
## 
## (1) demonstrate use of the algorithm for a given D, N
## (2) perform asymptotic analysis for a grid of both D, N
## 


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

# ------------------------------------------------------------------------------

# see new_features.R file for thorough step through of each of the values
# computed in the approximation (still in progress)

# x11() # move graphics to separate window

# model settings ---------------------------------------------------------------

D  = 2
Omega = diag(1, D)
# N  = 500
Sigma = diag(1, D)
Sigma_inv = solve(Omega)
mu_0 = rep(0, D)
prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, mu_0 = mu_0)

# closed form of the normalizing constant
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) # -3.683584, for D = 2, N = 500

# ------------------------------------------------------------------------------


# (1) demonstrate use of the algorithm for a given D, N ------------------------

set.seed(1)
J = 5000
N_approx = 1
u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
u_df_full = preprocess(u_samps, D, prior)
approx_skew = approx_lil(N_approx, D, u_df_full, J / N_approx, prior)
mean(approx_skew)

# ------------------------------------------------------------------------------

# (2) asymptotic analysis ------------------------------------------------------

N_vec_log = seq(6, 10, 0.02)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
logZ_0    = numeric(length(N_vec))  # store true value of log normalizing const
logZ      = numeric(length(N_vec))  # store approx of log normalizing const

print(length(N_vec))                # number of different sample sizes


for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    Sigma = D / N * Omega 
    Sigma_inv = solve(Sigma)
    
    mu_0 = rep(0, D)
    
    prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, mu_0 = mu_0)
    
    
    # compute true log normalizing const for given D, N
    logZ_0[i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma)

    
    print(paste("iter = ", i, "/", length(N_vec), 
                " -- Calculating LIL for D = ", D, ", N = ", N, sep = ''))
    
    
    J = 1e4         # number of total MC samples
    N_approx = 10   # number of approximations to form
    
    
    u_samps = rmvnorm(J, mean = rep(0, D), sigma = Sigma) %>% data.frame 
    u_df_full = preprocess(u_samps, D, prior = prior)
    approx_mvn = approx_lil(N_approx, D, u_df_full, J/N_approx, prior)
    
    logZ[i] = mean(approx_mvn)
    
} # end of loop iterating over different sample sizes



# plot the results so we can compute the slope of logZ ~ logN ------------------

lil_df = data.frame(logZ_0 = logZ_0, logZ = logZ, logn = log(N_vec))

lil_df = lil_df[is.finite(lil_df$logZ),] # omit values that have overflowed

lil_df %>% dim

lil_df_long = melt(lil_df, id.vars = "logn")

formula1 = y ~ x

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point() + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")




