


q = 3
p = 2
r = 1
D = r * p + r * q
N = p

## model : Y = ab' + E
# Y is (p x q)
# a is (p x 1)
# b is (q x 1)
# E is (p x q)

a_0 = rnorm(p) # (p x 1)
b_0 = rnorm(q) # (q x 1)

c_0 = a_0 %*% t(b_0)

nMCMC = 500       # number of MCMC samples from the posterior AFTER burnin
nBurn = 600       # number of samples to discard


rank1_obj = sampleRRR(nMCMC, nBurn, a_0, b_0, p, q, r, r, D, N, sig2, del, T)


u_samps = rank1_obj$u_samps # (nMCMC x D)

u = u_samps[300,] %>% unname %>% unlist

a_post = u[1:p]      # (p x 1)
b_post = tail(u, q)  # (q x 1)

(c_post = a_post %*% t(b_post))

c_0








