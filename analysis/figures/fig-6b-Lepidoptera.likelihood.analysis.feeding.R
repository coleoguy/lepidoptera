# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# read in data
load("../results/8b.calculate.likelihood.Lepidoptera.tip.feeding.RData")
liks <- unlist(x)
# fill in likelihood table
lik.table$likelihood <- liks
# plot the results
plot(liks,
     pch = 16,
     cex = .5,
     xlab = "Dropped tip",
     ylab = "Likelihood",
     col = c("gray", "red")[(liks > quantile(x = liks, 0.95)) + 1])
# cutoff mark
abline(h = quantile(x = liks, 0.95),
       col = "black",
       lty = 2,)
text(x = 420,
     y = quantile(x = liks, 0.95)-0.4,
     label = expression(paste("95"^"th","Quantile","")),
     pos = 3,
     cex = .7)
abline(h = full.tree.lik,
       col = "gray",
       lty = 2)
text(x = 420,
     y = full.tree.lik - .4,
     label = "Likelihood of the unpruned tree",
     pos = 3,
     cex = .7)
