
# ------------------------------------------------------------------------------


library(mvtnorm)           # for draws from multivariate normal
library("numDeriv")        # for grad() function - numerical differentiation
library('MCMCpack')        # for rinvgamma() function
library(sn)                # for rmns() function
library(VGAM)              # 
library(reshape2)          # for melt() function


# DELL_PATH = "C:/Users/chuu/mlike_approx"
LEN_PATH  = "C:/Users/ericc/mlike_approx"
# path for lenovo
setwd(LEN_PATH)

# path for dell
# setwd(DELL_PATH)

source("partition/partition.R")         # load partition extraction functions
source("hybrid_approx.R")               # load main algorithm functions
source("skew/mv_skew_normal_helper.R")  # load psi(), lambda()
source("extractPartition.R")


# fixed settings ---------------------------------------------------------------
D = 4
alpha = rep(1, D) 
mu_0 = rep(0, D)
Omega = diag(1, D)

# Sigma = D / N * Omega 
# Sigma_inv = solve(Sigma)

N_vec_log = seq(6, 10, 0.02)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data
logZ_0    = numeric(length(N_vec))  # store true value of log normalizing const

print(length(N_vec))                # number of different sample sizes


approx_taylor = matrix(NA, N_approx, length(N_vec))
approx_hybrid = matrix(NA, N_approx, length(N_vec))
approx_const  = matrix(NA, N_approx, length(N_vec))

for (i in 1:length(N_vec)) {
    
    N = N_vec[i]
    Sigma = D / N * Omega 
    Sigma_inv = solve(Sigma)
    
    # alpha, mu_0 initialized outside of loop, fixed for all values of N
    prior = list(Sigma = Sigma, Sigma_inv = Sigma_inv, 
                 alpha = alpha, mu_0 = mu_0)
    
    # compute true log normalizing const for given D, N
    logZ_0[i] = D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5)
    
    print(paste("iter = ", i, "/", length(N_vec), 
                " -- Calculating LIL for D = ", D, ", N = ", N, sep = ''))
    
    
    J = 1e4         # number of total MC samples
    N_approx = 50   # number of approximations to form
    u_samps = rmsn(J, xi = mu_0, Omega = Sigma, alpha = alpha) %>% data.frame 
    u_df_full = preprocess(u_samps, D, prior = prior)
    
    # approx_skew = approx_lil(N_approx, D, u_df_full, J/N_approx, prior)
    
    hml_skew = hml(N_approx, D, u_df_full, J / N_approx, prior)
    
    approx_taylor[,i] = hml_skew$taylor_vec
    approx_hybrid[,i] = hml_skew$hybrid_vec
    approx_const[,i]  = hml_skew$const_vec
    
    
} # end of loop iterating over different sample sizes


# true LIL ---------------------------------------------------------------------
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



# approx LIL -------------------------------------------------------------------

library(reshape2)

logZ_taylor = colMeans(approx_taylor)  # 
logZ_hybrid = colMeans(approx_hybrid)  #
logZ_const  = colMeans(approx_const)   #
logn        = log(N_vec)

lil_df = data.frame(logZ_0 = logZ_0, logZ_taylor = logZ_taylor, 
                    logZ_hybrid = logZ_hybrid, logZ_const = logZ_const, 
                    logn = logn)

lil_df_long = melt(lil_df, id.vars = "logn")


formula1 = y ~ x

x11()

ggplot(lil_df_long, aes(x = logn, y = value, 
                        color = as.factor(variable))) + geom_point(size = 0.7) + 
    geom_smooth(method = lm, se = F, formula = formula1) +
    labs(x = "log(n)", y = "log(Z)", 
         title = "Skew Normal (D = 2), True (Red), Hybrid (Blue),
         Taylor (Green), Constant (Purple)") + 
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none")




# overlays of both logZ, approx logZ -------------------------------------------

lil_df = data.frame(logZ_0 = logZ_0, logZ = logZ, logn = log(N_vec))

lil_df = lil_df[is.finite(lil_df$logZ),] # omit values that have overflowed

lil_df %>% dim

lil_df_long = melt(lil_df, id.vars = "logn")

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

rbind(logZ_0, logZ)[,1:7]










