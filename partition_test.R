
library(rpart)

gmix2_tree = rpart(psi_u ~ u1 + u2, data = u_df)

plotmo(gmix2_tree)

plotmo(gmix2_tree, nresponse = "present")

plot(gmix2_tree)
text(gmix2_tree, pretty = 0, cex = 0.7)


par(mfrow = c(1,2))

plot(u_tree)
text(u_tree, pretty = 0, cex = 0.7)

# overlay partition on scatterplot of points drawn from true density
plot(u_df[,1], u_df[,2], pch = 20, cex = 0.8, col = "pink",
     xlab = 'u1', ylab = 'u2', main = 'pi1 = 0.2, pi2 = 0.8 J = 5000')
partition.tree(u_tree, add = TRUE, cex = 0.0001, ordvars = c("u1", "u2"))


par(mfrow = c(1,1))
