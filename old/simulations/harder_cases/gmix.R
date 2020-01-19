

## preliminary steps to check if the partitioning is doing the right thing

# -----------------------------------------------------------------------------

library('mvtnorm')   # multivariate normal density
library('MASS')      # mvnorm()
library('ggplot2')
library('tree')

set.seed(123)

d = 2
N = 100
mu1 = c(1, 2)
mu2 = c(4, 6)
Sigma1 = d / N * diag(c(1, 1))
Sigma2 = d / N * diag(c(1, 1))

# mixture weights
pi1 = 0.2
pi2 = 0.8

J = 5000

n1 = sum(runif(J) <= pi1) # draws from 1st gaussian
n2 = J - n1               # draw sfrom 2nd gaussian


# define psi() function --------------------------------------------------------
psi = function(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2) {
    
    # log-density evaluated over the rows of u -- (J x 1)
    loglik = log(pi1 * dmvnorm(u, mu1, Sigma1) + 
                     pi2 * dmvnorm(u, mu2, Sigma2))
    
    return(loglik)
    
} # end psi() function


u = rbind(mvrnorm(n1, mu1, Sigma1), mvrnorm(n2, mu2, Sigma2)) # (J x 2)

# check distribution
u_df = data.frame(u1 = u[,1], u2 = u[,2])
# ggplot(u_df, aes(u1, u2)) + geom_point() + theme_bw()

# evaluate psi(u)
psi_u = psi(u, mu1, mu2, Sigma1, Sigma2, pi1, pi2)           # (J x 1)
u_df = data.frame(u1 = u_df$u1, u2 = u_df$u2, psi_u = psi_u) # (J x 3)

# fit decision tree
u_tree = tree(psi_u ~ u1 + u2, u_df)

# overlay partition on scatterplot of points drawn from true density
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.8, col = "pink",
     xlab = 'u1', ylab = 'u2', main = 'pi1 = 0.2, pi2 = 0.8 J = 5000')
partition.tree(u_tree, add = TRUE, cex = 0.0001, ordvars = c("u1", "u2"))





