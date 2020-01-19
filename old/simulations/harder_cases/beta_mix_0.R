
library('ggplot2')
library('tree')

# beta parameters
alpha1 = 2
beta1 = 10
alpha2 = 8
beta2 = 2

hist(rbeta(1000, 2, 10))
hist(rbeta(1000, 8, 2))

set.seed(123)

# mixture weights
pi1 = 0.2
pi2 = 0.8

J = 10000

n1 = sum(runif(J) <= pi1) # draws from 1st gaussian
n2 = J - n1

# define psi() function --------------------------------------------------------
psi = function(u, alpha1, alpha2, beta1, beta2) {
    
    loglik = log(pi1 * dbeta(u, alpha1, beta1) + pi2 * dbeta(u, alpha2, beta2))
    
    return(loglik)    
    
} # end psi() function

# sample from the mixture density
u1 = rbeta(n1, alpha1, beta1)
u2 = rbeta(n2, alpha2, beta2)
u = c(u1, u2)   # (J x 1)
hist(c(u1, u2)) # mixture of beta distribution

psi_u = psi(u, alpha1, alpha2, beta1, beta2)
u_df = data.frame(u = u, psi_u = psi_u) # (J x 3)

# make minsize smaller to encourage a larger tree
u_tree0 = tree(psi_u ~ u, u_df)
u_tree1 = tree(psi_u ~ u, u_df, control = tree.control(nobs = J, mincut = 1, 
                                                       minsize = 2))

# plot the tree with labeled branches
plot(u_tree0)
text(u_tree0, all = TRUE, cex = 0.7)
plot(u_tree1)
text(u_tree1, all = TRUE, cex = 0.7)

par(mfrow = c(1,2))
# predicted value is plotted against the predictor over the range in the 
# training set
partition.tree(u_tree0, add = F, cex = 1.5)
partition.tree(u_tree1, add = F, cex = 1.5)
hist(c(u1, u2)) # mixture of beta distribution

par(mfrow = c(1,1))







