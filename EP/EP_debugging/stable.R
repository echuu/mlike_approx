

lpnorm = function(x) {
  pnorm(x, log.p = TRUE)
}

Z_i_hat = exp(lpnorm(a * sqrt(2))) * 
  expm1(lpnorm(b * sqrt(2)) - lpnorm(a * sqrt(2)))

mu_i_hat = mu + sqrt(sigma / (2 * pi)) * exp(-b^2 - lpnorm(a * sqrt(2))) * 
  expm1(-a^2 + b^2) / expm1(lpnorm(b * sqrt(2)) - lpnorm(a * sqrt(2)))

mu_hat

sigma_i_hat = mu^2 + sigma - mu_i_hat^2 + sqrt(sigma / (2 * pi)) * 
  1 / Z_i_hat * ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))

sigma_i_hat

sigma_hat

all.equal(sigma_i_hat, Re(sigma_hat))
sigma_i_hat == Re(sigma_hat)


sourceCpp("C:/Users/ericc/rcpp-epmgp/src/axisepmgp.cpp")

res_pkg = axisepmgp(m, K, lb, ub)

res_pkg$logZ
res_pkg$mu
res_pkg$Sigma

tn_moments = function(lb_vec, ub_vec, mu_in, sigma_in) {

  d = length(lb_vec)
  logz_hat_out  = numeric(d)
  z_hat_out     = numeric(d)
  mu_hat_out    = numeric(d)
  sigma_hat_out = numeric(d)
  
  for (i in 1:d) {
    
    lb    = lb_vec[i]
    ub    = ub_vec[i]
    mu    = mu_in[i]
    sigma = sigma_in[i]
    
    logz_hat_other_tail = 0
    logz_hat = 0
    mean_const = 0
    var_const = 0
    
    ## establish bounds
    a = (lb - mu) / sqrt(2 * sigma)
    b = (ub - mu) / sqrt(2 * sigma)
    
    Z_i_hat = exp(lpnorm(a * sqrt(2))) * 
      expm1(lpnorm(b * sqrt(2)) - lpnorm(a * sqrt(2)))
    
    mu_i_hat = mu + sqrt(sigma / (2 * pi)) * exp(-b^2 - lpnorm(a * sqrt(2))) * 
      expm1(-a^2 + b^2) / expm1(lpnorm(b * sqrt(2)) - lpnorm(a * sqrt(2)))
    
    sigma_i_hat = mu^2 + sigma - mu_i_hat^2 + sqrt(sigma / (2 * pi)) * 
      1 / Z_i_hat * ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))
    
    z_hat     = Z_i_hat
    
    logz_hat_out[i]  = log(z_hat)  
    z_hat_out[i]     = z_hat   
    mu_hat_out[i]    = mu_i_hat  
    sigma_hat_out[i] = sigma_i_hat 
    
  } # end of for loop
  
  result = list(logz_hat  = (logz_hat_out), 
                z_hat     = (z_hat_out), 
                mu_hat    = (mu_hat_out), 
                sigma_hat = (sigma_hat_out))
  return(result)
}



tn_cpp = trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity)
tn_r   = tn_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity)


all.equal(
  lapply(tn_cpp, as.matrix),
  lapply(tn_r, as.matrix)
)




