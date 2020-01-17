


## TODO: try the outlined code below using the MVN-IG example

## we should be able to still pre-process everything -- pre-determine the
## K_sims samples (process should still be very similar)

# outline for analyzing LIL_N vs. log(N)

LIL_N = numeric(length)

# values of N for which we will compute + approximate the LIL
N_vec  = c(50, 100, 200, 500, 750, 1000, 2000, 4000, 8000, 10000)

K_sims      = 200 # number of simulations to run FOR EACH N in N_vec
M_estimates = 1   # num of (hybrid) estimates to calc FOR EACH of the K_sims


for (i in 1:length(N_vec)) {
    
    N = N[i]
    
    LIL_N_k = numeric(K_sims) # store the K_sims results
    
    for (k in 1:K_sims) {
        
        ## simulate N data points + sample from posterior
        
        ## form the hybrid approximation
        
        
    }
    
    LIL_N[i] = mean(LIL_N_k)
    
}


LIL_df = data.frame(LIL_N = LIL_N, log_N = log(N_vec))
plot(LIL_N ~ log_N)

lm(LIL_N ~ log_N, LIL_df)









