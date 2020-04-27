





# TODO: for each partition, find u that corresponds to max value of psi




# for each partition, find u that most closely corresponds to psi_star
psi_df = hml_approx$u_df_fit %>% mutate(psi_dev = abs(psi_u - psi_star))

psi_star_df = psi_df %>% 
    group_by(leaf_id) %>% 
    slice(which.min(psi_dev)) %>%  # extract rows that minimize log(Q(c))
    data.frame()

x11()
plotPartition(u_df, hml_approx$param_out, psi_df, psi_star_df)

## add text for psi_max value
## add text for psi labels for each point (smaller than the psi_star labels)


plotPartition = function(u_df, param_out, psi_df, psi_star_df, OVERLAY = TRUE) {
    
    plot(u_df[,1], u_df[,2], pch = 20, cex = 1, col = "white",
         xlab = 'u1', ylab = 'u2', main = '')
    
    rect(param_out$u1_lb, 
         param_out$u2_lb, 
         param_out$u1_ub, 
         param_out$u2_ub)
    
    # add psi_hat labels for each partition
    # text(x = param_out$u1_lb + (param_out$u1_ub - param_out$u1_lb) / 2, 
    #      y = param_out$u2_lb + (param_out$u2_ub - param_out$u2_lb) / 2,
    #      labels = round(param_out$psi_star, 5),
    #      col = 'purple',
    #      cex = 1.5)
    
    # add psi values for all points on the plot
    text(x = psi_df$u1, 
         y = psi_df$u2,
         labels = round(psi_df$psi_u, 3),
         cex = 0.5)
    
    text(x = psi_star_df$u1, y = psi_star_df$u2,
         labels = round(param_out$psi_star, 3),
         col = 'white',
         cex = 0.5)
    
    text(x = psi_star_df$u1, y = psi_star_df$u2,
         labels = round(param_out$psi_star, 3),
         col = 'purple',
         cex = 1.2)
    
    # make the 'median' points red and large
    if (OVERLAY) {
        points(x = psi_star_df$u1, y = psi_star_df$u2,
               col = 'red', pch = 19, cex = 1.2)
    }
    
}


plotPartition = function(u_df, param_out, psi_df, psi_star_df, OVERLAY = FALSE) {
    
    plot(u_df[,1], u_df[,2], pch = 20, cex = 1, col = "cyan",
         xlab = 'u1', ylab = 'u2', main = '')
    
    rect(param_out$u1_lb, 
         param_out$u2_lb, 
         param_out$u1_ub, 
         param_out$u2_ub)
    
    # add psi_hat labels for each partition
    # text(x = param_out$u1_lb + (param_out$u1_ub - param_out$u1_lb) / 2, 
    #      y = param_out$u2_lb + (param_out$u2_ub - param_out$u2_lb) / 2,
    #      labels = round(param_out$psi_star, 5),
    #      col = 'purple',
    #      cex = 1.5)
    
    # text(x = psi_star_df$u1, y = psi_star_df$u2,
    #      labels = round(param_out$psi_star, 5),
    #      col = 'purple',
    #      cex = 1.5)
    
    # add psi values for all points on the plot
    # text(x = psi_df$u1, 
    #      y = psi_df$u2,
    #      labels = round(psi_df$psi_u, 5),
    #      cex = 0.5)
    
    # make the 'median' points red and large
    if (OVERLAY) {
        points(x = psi_star_df$u1, y = psi_star_df$u2,
               col = 'red', pch = 19, cex = 1.2)
    }
    
}




