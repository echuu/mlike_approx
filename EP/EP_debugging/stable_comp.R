

lpnorm = function(x) {
  pnorm(x, log.p = TRUE)
}

## run these to manually set values of a, b
a = 3
b = 5



## run these to get original value of a, b
a = (lb - mu) / sqrt(2 * sigma)  
b = (ub - mu) / sqrt(2 * sigma)


b = b + 1e-9
a
b


a = -0.005
b = -0.004


all.equal(a, b)

exp_a2b2 = exp(a^2 - b^2)
exp_b2a2 = exp(b^2 - a^2)
erfcx_a = RcppFaddeeva::erfcx(abs(a));
erfcx_b = RcppFaddeeva::erfcx(abs(b));


2 * exp(a^2 - b^2 + pnorm(-b * sqrt(2), log.p = T)) * 
  expm1(pnorm(-a * sqrt(2), log.p = T) - pnorm(-b * sqrt(2), log.p = T) + b^2)

2 * exp(a^2 - b^2 + pnorm(-b * sqrt(2), log.p = T)) * 
  (exp(pnorm(-a * sqrt(2), log.p = T) - pnorm(-b * sqrt(2), log.p = T) + b^2 ) - 1)


#### typical case 3 -- sign(a) == sign(b): 

## mean_const calculations ---------------

erfcx_a - exp_a2b2 * erfcx_b ## original

(x1 = 2 * exp(a^2 + lpnorm(- abs(b) * sqrt(2))) * 
    expm1(lpnorm(-abs(a) * sqrt(2)) - lpnorm(-abs(b) * sqrt(2))))

exp_b2a2 * erfcx_a - erfcx_b ## original

(x2 = 2 * exp(b^2 + lpnorm(- abs(b) * sqrt(2))) * 
    expm1(lpnorm(-abs(a) * sqrt(2)) - lpnorm(-abs(b) * sqrt(2))))


erfcx_a
exp(a^2) * 2 * pnorm(-abs(a) * sqrt(2))
erfcx_b
exp(b^2) * 2 * pnorm(-abs(b) * sqrt(2))


## entire calculation:
2 * sign(a) * (
  1 / ( erfcx_a - exp_a2b2 * erfcx_b ) -
    1 / ( exp_b2a2 * erfcx_a - erfcx_b )
)

## 'stable' calculation:
mean_const_stable = 2 * sign(a) * (1/x1 - 1/x2)
mean_const_stable
## var_const calculations ---------------


## var_const: entire original calculation: 
2 * sign(a) * (
  (lb + mu) / (erfcx_a - exp_a2b2 * erfcx_b) - 
  (ub + mu) / (exp_b2a2 * erfcx_a - erfcx_b)
)

## 'stable' calculation:
var_const_stable = 2 * sign(a) * ((lb + mu) / x1 - (ub + mu) / x2)
var_const_stable


#### parameter updates:
mu_hat_stable = mu + mean_const_stable * sqrt(sigma / (2 * pi))
sigma_hat_stable = sigma + mu^2 - mu_hat_stable^2 + 
  var_const_stable * sqrt(sigma / 2 * pi)
sigma_hat_stable


z = exp(pnorm(a * sqrt(2), log.p = TRUE)) * 
  expm1(pnorm(b * sqrt(2), log.p = TRUE) - pnorm(a * sqrt(2), log.p = TRUE))
z

other_var = mu^2 + sigma - mu_hat_stable^2 + 1/z * sqrt(sigma / (2 * pi)) * 
  ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))
other_var

# 1/z *  ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))

1/other_var - tau_cavity - tau_site

mu_hat
sigma_hat

mu_hat_stable
sigma_hat_stable










