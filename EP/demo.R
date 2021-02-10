

library(BDgraph)

# Eq. (2.4) in Uhler paper (p.12)
I_G = function(delta) {
    7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) + 
        lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) + 
        lgamma(delta + 5/2)
}


d = 5
delta = 3
# |E| = 7. Applying formula at the end of pg.2 
# of the Uhler paper. 
#log normalizing  constant of the G-WIshart is 
log(2^(0.5*d*delta + 7)) + I_G(0.5*(delta-2))


d = 5 # number vertices
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), d, d)


I_G(delta)
V = diag(1, d)
gnorm(G_5, b, V, 100)

