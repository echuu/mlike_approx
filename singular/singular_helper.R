


# psi() : negative log posterior
psi = function(u, N) {
    
    return(N * u[1]^2 * u[2]^4)
    
}


# lambda() : gradient of psi
lambda = function(u, N) {
    
    l_1 = 2 * N * u[1] * u[2]^4
    l_2 = 4 * N * u[1]^2 * u[2]^3
    
    return(c(l_1, l_2))
}


preprocess = function(stan_fit, D, N) {
    
    u_samp = rstan::extract(gamma_fit, pars = c("u"), permuted = TRUE)
    
    u_post = u_samp$u %>% data.frame() # (J * N_approx) x 2
    
    
    psi_u = apply(u_post, 1, psi, N = N) %>% unname() # (J x 1)
    
    # (1.2) construct u_df -- this will require some automation for colnames
    u_df_names = character(D + 1)
    for (d in 1:D) {
        u_df_names[d] = paste("u", d, sep = '')
    }
    u_df_names[D + 1] = "psi_u"
    
    # populate u_df
    u_df = cbind(u_post, psi_u) # (J * N_approx) x (D + 1)
    
    # rename columns (needed since these are referenced explicitly in partition.R)
    names(u_df) = u_df_names
    
    
    return(u_df)
    
}








