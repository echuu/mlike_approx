


# compute normalizing constant for G_5

# Eq. (2.4) in Uhler paper (p.12)
I_G = function(delta) {
    7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) + 
        lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) + 
        lgamma(delta + 5/2)
}

delta = b
I_G(delta)

p = 5 # number vertices
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)


# target quantity is Eq. (2.2) on p.9 of Uhler




