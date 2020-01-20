
library(ggplot2)

set.seed(1)

mu = 30
sigmasq = 4
m_0 = 0 
tau = 3


N_grid = c(50, 100, 200, 500, 1000, 10000)
iters = length(N_grid)


lil = numeric(iters)
lil_n = numeric(1000)

set.seed(1)
for (i in 1:iters) {
    
    
    N = N_grid[i]
    
    # compute log marginal likelihood 1000 times
    for (t in 1:length(lil_n)) {
        
        x = rnorm(N, mean = mu, sd = sqrt(sigmasq))
        xbar = mean(x)
        
        lil_n[t] = log(sqrt(sigmasq)) - N * log(2 * pi * sqrt(sigmasq)) - 
            0.5 * log(N * tau^2 + sigmasq) - sum(x^2) / (2 * sigmasq) - 
            m_0^2 / (2 * tau^2) + 1 / (2 * (N * tau^2 + sigmasq)) * 
            (tau^2 * N^2 * xbar^2 / sigmasq + sigmasq * m_0^2 / tau^2 + 
                 2 * N * xbar * m_0)
        
    }
    
    lil[i] = mean(lil_n)
    
}

lil_df = data.frame(n = N_grid, lil = lil)
ggplot(lil_df, aes(log(n), lil)) + geom_point() + theme_bw()


ggplot(lil_df, aes(n, lil)) + geom_point() + theme_bw()
lm(lil ~ n, lil_df)







