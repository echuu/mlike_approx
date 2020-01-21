

install.packages("pracma")
library(pracma)
fun <- function(x, y) exp(-n*x^2*y^4)
logn = seq(1,20,by=0.1)
N = exp(logn)
lML = rep(0,length(N))

i = 1;
n = 100
n = 1000
for (n in N){
    result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
    lML[i] = log(result$Q)
    i = i + 1;
}
plot(logn,lML)
lm(lML ~ logn)


n = 100
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -0.724

n = 1000
result = integral2(fun, 0, 1, 0, 1, reltol = 1e-50)
log(result$Q) # -1.223014


N = 1000
D = 2                               # dimension of parameter

gamma_dat = list(N = N)

# should give us (J * N_approx) draws
gamma_fit_N = stan(file    =  'gamma_sample.stan', 
                   data    =  gamma_dat,
                   iter    =  J_iter,
                   warmup  =  burn_in,
                   chains  =  n_chains,                  
                   seed    =  stan_seed,
                   control =  list(adapt_delta = 0.99),  
                   refresh = 0)          

u_df_N = preprocess(gamma_fit_N, D, N)

rpart_sing = rpart(psi_u ~ ., u_df_N, cp = 0.001)
plot(rpart_sing)
text(rpart_sing, size = 0.5)

approx_test = approx_lil_stan(N_approx, D, N, u_df_N, J)

approx_test

mean(approx_test)

# ------------------------------------------------------------------------------

u_df_N %>% head

u_rpart = rpart(psi_u ~ ., u_df_in)

rpart_obj = u_rpart


u_rpart = sing_rpart

 
plot(u_df_N[,1], u_df_N[,2])

plot(sing_rpart)
text(sing_rpart, size = 0.5)




D = 10
N = 7 # pseudo-sample size
Omega = diag(1, D)
Sigma = D / N * Omega 
Sigma_inv = solve(Sigma)
alpha = rep(1, D) 
mu_0 = rep(0, D)

# -4.3767 for D = 2, N = 500
D / 2 * log(2 * pi) + 0.5 * log_det(Sigma) + log(0.5) 









