



library('mvtnorm')    # multivariate normal density
library('MASS')       # mvnorm()
library('ggplot2')
library('tree')
library('data.table') # data management

set.seed(123)

d = 2

mu1 = c(1, 2)
mu2 = c(4, 6)
mu3 = c(4, 2)
Sigma1 = d / N * diag(c(1, 1))
Sigma2 = d / N * diag(c(1, 1))
Sigma3 = d / N * diag(c(1, 1))

pi1 = 0.2
pi2 = 0.5
pi3 = 0.3

J = 5000
draws = runif(J)

n1 = sum((draws <= pi1))
n2 = sum((draws > pi1) & (draws <= pi1 + pi2))
n3 = sum((draws > (pi1 + pi2)))

# define psi() function --------------------------------------------------------
psi = function(u, mu1, mu2, Sigma1, Sigma2, Sigma3, pi1, pi2, pi3) {
    
    # log-density evaluated over the rows of u -- (J x 1)
    loglik = log(pi1 * dmvnorm(u, mu1, Sigma1) + 
                 pi2 * dmvnorm(u, mu2, Sigma2) + 
                 pi3 * dmvnorm(u, mu3, Sigma3))
    
    return(loglik)
    
} # end psi() function

u_list = list(mvrnorm(n1, mu1, Sigma1), 
              mvrnorm(n2, mu2, Sigma2),
              mvrnorm(n3, mu3, Sigma3)) 

u = rbindlist(lapply(u_list, as.data.frame))
u_df = data.frame(u1 = u[,1], u2 = u[,2])
names(u_df) = c("u1", "u2")
# ggplot(u_df, aes(u1, u2)) + geom_point() + theme_bw()


# evaluate psi(u)
psi_u = psi(u, mu1, mu2, Sigma1, Sigma2, Sigma3, pi1, pi2, pi3)   # (J x 1)
u_df = data.frame(u1 = u_df$u1, u2 = u_df$u2, psi_u = psi_u)      # (J x 3)

# fit decision tree
u_tree = tree(psi_u ~ u1 + u2, u_df)

# overlay partition on scatterplot of points drawn from true density
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.8, col = "green",
     xlab = 'u1', ylab = 'u2', 
     main = 'pi1 = 0.2, pi2 = 0.5, pi3 = 0.3, J = 1000')
partition.tree(u_tree, add = TRUE, cex = 0.0001, ordvars = c("u1", "u2"))



