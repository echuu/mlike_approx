


sigma + var_const * sqrt(sigma / (2 * pi)) + mu^2 - mu_hat^2  

sign(a)
sign(b)

abs(a)
abs(b)



# mu_hat update
mu + 1/tmp1 * sqrt(sigma / (2*pi)) * (exp(-a^2) - exp(-b^2))

# mean_const value
1/tmp1 * (exp(-a^2) - exp(-b^2))


## sigma_hat update
mu^2 - mu_hat^2 + sigma + 1 / tmp1 * sqrt(sigma / (2 * pi)) * 
  ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))

sigma + var_const * sqrt(sigma / (2 * pi)) + mu^2 - mu_hat^2  

1/tmp1 * ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))

b > a
sign(a)




2 * sign(a) * (
  1 / (( erfcx_a - exp_a2b2 * erfcx_b )) -
    1 / (( exp_b2a2 * erfcx_a - erfcx_b )))

erfcx_a - exp_a2b2 * erfcx_b

exp(a^2 + pnorm(-sqrt(2) * abs(b), log.p = T)) * 
  expm1(pnorm(-sqrt(2) * abs(a), log.p = TRUE) - pnorm(-sqrt(2) * abs(b), log.p = TRUE))

# ------------------------------------------------------------------------------
a = -5
b = -3

exp_a2b2 = exp(a^2 - b^2)
exp_b2a2 = exp(b^2 - a^2)
erfcx_a = RcppFaddeeva::erfcx(abs(a));
erfcx_b = RcppFaddeeva::erfcx(abs(b));

mean_const = 2 * sign(a) * (
  1 / (( erfcx_a - exp_a2b2 * erfcx_b )) -
    1 / (( exp_b2a2 * erfcx_a - erfcx_b )))
mean_const


z = exp(pnorm(a * sqrt(2), log.p = TRUE)) * 
  expm1(pnorm(b * sqrt(2), log.p = TRUE) - pnorm(a * sqrt(2), log.p = TRUE))
z1 = 1/2 * (pchisq(2*b^2, 1) * sign(b) - pchisq(2*a^2, 1) * sign(a))
z2 = 0.5 * exp(pchisq(2*b^2, 1, log.p = TRUE)) * 
  expm1(pchisq(2*a^2, 1, log.p = TRUE) - pchisq(2*b^2, 1, log.p = TRUE))

z
z1

# compare this to mu_hat in orig code
(mu_hat_0 = mu + 1 / z * sqrt(sigma / (2 * pi)) * exp(-b^2) * expm1(-a^2 + b^2))

1 / z * exp(-b^2) * expm1(-a^2 + b^2)

mean_const

mu + mean_const * sqrt(sigma / (2 * pi))
(mu + sqrt(sigma / (2 * pi)) * 1 / z * exp(-b^2) * expm1(-a^2 + b^2))^2


# ------------------------------------------------------------------------------


lpnorm = function(x) {
  pnorm(x, log.p = TRUE)
}

### expressions for all 3 updates:

# Z_i 
z = exp(pnorm(a * sqrt(2), log.p = TRUE)) * 
  expm1(pnorm(b * sqrt(2), log.p = TRUE) - pnorm(a * sqrt(2), log.p = TRUE))

mu_i_hat = mu + sqrt(sigma / (2 * pi)) / z * exp(-b^2) * expm1(-a^2 + b^2)

sigmasq_i_hat = mu^2 - mu_i_hat^2 + sigma + sqrt(sigma / (2 * pi)) / z * 
  ((lb + mu) * exp(-a^2) - (mu + ub) * exp(-b^2))

mu_hat
mu_i_hat

sigma_hat
sigmasq_i_hat


a = 3
b = 5
a = (lb - mu) / sqrt(2 * sigma)
b = (ub - mu) / sqrt(2 * sigma)
exp_a2b2 = exp(a^2 - b^2)
exp_b2a2 = exp(b^2 - a^2)
erfcx_a = RcppFaddeeva::erfcx(abs(a));
erfcx_b = RcppFaddeeva::erfcx(abs(b));




2 * exp(a^2 - b^2 + pnorm(-b * sqrt(2), log.p = T)) * 
  expm1(pnorm(-a * sqrt(2), log.p = T) - pnorm(-b * sqrt(2), log.p = T) + b^2)

2 * exp(a^2 - b^2 + pnorm(-b * sqrt(2), log.p = T)) * 
  (exp(pnorm(-a * sqrt(2), log.p = T) - pnorm(-b * sqrt(2), log.p = T) + b^2 ) - 1)


## typical case 3 -- sign(a) == sign(b): mean_const calculations ---------------

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
2 * sign(a) * (1/x1 - 1/x2)



# ------------------------------------------------------------------------------
a = -3
b = 4

a = (lb - mu) / sqrt(2 * sigma)
b = (ub - mu) / sqrt(2 * sigma)

exp_a2b2 = exp(a^2 - b^2)
exp_b2a2 = exp(b^2 - a^2)
erfcx_a = RcppFaddeeva::erfcx(abs(a));
erfcx_b = RcppFaddeeva::erfcx(abs(b));

z1 = 1/2 * (pchisq(2*b^2, 1) * sign(b) - pchisq(2*a^2, 1) * sign(a))
z2 = 0.5 * exp(pchisq(2*b^2, 1, log.p = TRUE)) * 
  expm1(pchisq(2*a^2, 1, log.p = TRUE) - pchisq(2*b^2, 1, log.p = TRUE))
z3 = 
z1 == z2


mu + sqrt(sigma / (2* pi))

# mean_const
2 * sign(a) * 
  (1 / ( erfcx_a - exp_a2b2 * erfcx_b ) - 1 / ( exp_b2a2 * erfcx_a - erfcx_b ))

# mean_const as derived from the paper:
1 / z2 * (exp(-a^2) - exp(-b^2))
1 / z2 * (exp(-b^2) * expm1(-a^2+b^2))

(exp(-b^2) * expm1(-a^2+b^2)) == ((exp(-a^2) - exp(-b^2)))

## mu_hat given in code
mu + mean_const * sqrt(sigma / (2 * pi))

### expression given in the paper for mu_hat update
mu_hat_0 = mu + 1 / tmp1 * sqrt(sigma / (2 * pi)) * (exp(-a^2) - exp(-b^2))
mu_hat_0






# var_const as derived from the paper:
mu^2 - mu_hat_0^2 + sigma + 1 / tmp1 * sqrt(sigma / (2 * pi)) * 
  ((lb + mu) * exp(-a^2) - (ub + mu) * exp(-b^2))


# var_const given in the code (unstable version)
sigma + var_const * sqrt(sigma / (2 * pi)) + mu^2 - mu_hat_0^2  






