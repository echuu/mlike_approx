
set.seed(1)
mu = 30
sigmasq = 4
m_0 = 0 
tau = 3

N_grid = c(50, 100, 200, 500, 1000, 10000)
iters = length(N_grid)


n_sims = 200 # number of times to compute LIL for EACH value of N

lil_N = numeric(iters) # store LIL for each N in the grid

for (i in 1:iters) {
    
    N = N_grid[i]
    
    lil = numeric(n_sims) # store log marginal likelihood for each simulation
    
    for (j in 1:n_sims) {
        
        x = rnorm(N, mean = mu, sd = sqrt(sigmasq))
        xbar = mean(x)
        
        lil[j] = log(sqrt(sigmasq)) - N * log(2 * pi * sqrt(sigmasq)) - 
            0.5 * log(N * tau^2 + sigmasq) - sum(x^2) / (2 * sigmasq) - 
            m_0^2 / (2 * tau^2) + 1 / (2 * (N * tau^2 + sigmasq)) * 
            (tau^2 * N^2 * xbar^2 / sigmasq + sigmasq * m_0^2 / tau^2 + 
                 2 * N * xbar * m_0) - 
            sum(dnorm(x, mean = xbar, sd = sqrt(sigmasq), log = TRUE))
        
    }
    
    lil_N[i] = mean(lil)
    
}

lil_df = data.frame(log_n = log(N_grid), lil_n = lil_N)


ggplot(lil_df, aes(log_n, lil_n)) + geom_point()

lm(lil_n ~ log_n, lil_df)












