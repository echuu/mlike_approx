

source("C:/Users/ericc/mlike_approx/partition/partition.R")
setwd("C:/Users/ericc/mlike_approx/singular")
source("singular_helper.R")

library("ggplot2")
library("ggpmisc")



sim1 = read.csv("approx_N50_J1000_grid61.csv")

N_vec_log = seq(4, 10, 0.1)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data

log_z_n = colMeans(sim1)
log_n   = log(N_vec)

lil_df = data.frame(log_z_n, log_n)

formula1 = y ~ x


ggplot(lil_df, aes(log_n, log_z_n)) + geom_point() + 
    labs(title = "50 approximations (1000 MC samples each) for each N") + 
    geom_smooth(method = lm, se = T, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16)

lm(log_z_n ~ log_n, lil_df) # slope should be -0.25



sim2 = read.csv("approx_N50_J2000_grid121.csv")

N_vec_log = seq(4, 10, 0.05)        # sample size that is uniform over log scale
N_vec     = floor(exp(N_vec_log))   # sample size to use to generate data

log_z_n = colMeans(sim2)
log_n   = log(N_vec)

lil_df = data.frame(log_z_n, log_n)

formula1 = y ~ x


ggplot(lil_df, aes(log_n, log_z_n)) + geom_point() + 
    labs(title = "50 approximations (2000 MC samples each) for each N") + 
    geom_smooth(method = lm, se = T, formula = formula1) +
    stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "top",
                 eq.with.lhs = "logZ~`=`~",
                 eq.x.rhs = "~logN",
                 formula = formula1, parse = TRUE, size = 8) +
    theme_bw(base_size = 16)

lm(log_z_n ~ log_n, lil_df) # slope should be -0.25






