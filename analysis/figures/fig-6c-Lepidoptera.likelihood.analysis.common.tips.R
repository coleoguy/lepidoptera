# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# read in data
load("../results/8a.calculate.likelihood.Lepidoptera.tip.RData")
liks <- unlist(x)
# fill in likelihood table
lik.table$likelihood <- liks
lik.table.1 <- lik.table
full.tree.lik.1 <- full.tree.lik
# clear results but keep the like table
rm(list = ls()[c(-13,-27)])
# read in data
load("../results/8b.calculate.likelihood.Lepidoptera.tip.feeding.RData")
liks <- unlist(x)
# fill in likelihood table
lik.table$likelihood <- liks
lik.table.2 <- lik.table
full.tree.lik.2 <- full.tree.lik
# get the common tips
commonTips <- lik.table.1$likelihood > quantile(x = lik.table.1$likelihood, 0.95) & lik.table.2$likelihood > quantile(x = lik.table.2$likelihood, 0.95)

par(mfcol = c(1,2))
# plot the results
plot(lik.table.1$likelihood,
     pch = 16,
     cex = .5,
     xlab = "Dropped tip",
     ylab = "Likelihood",
     col = c("gray", "red")[(lik.table.1$likelihood > quantile(x = lik.table.1$likelihood, 0.95)) + 1])
#plot title
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)

# highlight the 
points(x = which(commonTips),
       y = lik.table.1$likelihood[commonTips],
       pch = 21,
       cex = 2,
       col = "black")

# cutoff mark
abline(h = quantile(x = lik.table.1$likelihood, 0.95),
       col = "black",
       lty = 2,)
text(x = 420,
     y = quantile(x = lik.table.1$likelihood, 0.95)-0.4,
     label = expression(paste("95"^"th","Quantile","")),
     pos = 3,
     cex = .7)
abline(h = full.tree.lik.1,
       col = "gray",
       lty = 2)
text(x = 440,
     y = full.tree.lik.1 - .4,
     label = "Likelihood of the unpruned tree",
     pos = 2,
     cex = .7)
# plot table 2
# plot the results
plot(lik.table.2$likelihood,
     pch = 16,
     cex = .5,
     xlab = "Dropped tip",
     ylab = "Likelihood",
     col = c("gray", "red")[(lik.table.2$likelihood > quantile(x = lik.table.2$likelihood, 0.95)) + 1])
#plot title
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# highlight the 
points(x = which(commonTips),
       y = lik.table.2$likelihood[commonTips],
       pch = 21,
       cex = 2,
       col = "black")

# cutoff mark
abline(h = quantile(x = lik.table.2$likelihood, 0.95),
       col = "black",
       lty = 2,)
text(x = 420,
     y = quantile(x = lik.table.2$likelihood, 0.95)-0.4,
     label = expression(paste("95"^"th","Quantile","")),
     pos = 3,
     cex = .7)
abline(h = full.tree.lik.2,
       col = "gray",
       lty = 2)
text(x = 440,
     y = full.tree.lik.2 - .4,
     label = "Likelihood of the unpruned tree",
     pos = 2,
     cex = .7)
# save the results
write.csv(commonTips, file = "../results/potential.tips.impacting.model.csv",row.names = F)
