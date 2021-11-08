# Terrence Sylvester
# 12 November 2020
# pradakshanas@gmail.com

# read in genetic diversity data
dat <- read.csv("../data/genetic.diversity/genetic.diversity.csv", as.is = T)
# plot
plot(x = NULL,
     y = NULL,
     xlim = c(.5,2.5),
     ylim = c(0,max(dat$X._0D)),
     xlab = "Host type",
     ylab = "Genetic diversity at four fold degenerative sites",
     axes = F)
axis(side = 1, at = c(.5,2.5),labels = F)
axis(side = 2)
boxplot(dat$X._0D ~ dat$LHP.breadth, add = T,axes = F, col = NULL)
points(x =jitter(c(1,2)[as.factor(dat$LHP.breadth)]),
       y = dat$X._0D,
       cex = 1,
       pch = 16)
axis(side = 1, at = c(1,2), labels = c("Generalist", "Specialist"))
# perform statistical test
t.test <- t.test(dat$geneticDiversity ~ dat$type)
text(x = 0.5,
     y = 0.035,
     labels = paste("t test p-value ", round(t.test$p.value, 2), sep = ""),
     pos = 4)
# 0 fold degenerate sites
# make a new datase which have genetic diversity and feeding type
dat0D <- as.data.frame(matrix(data = NA,
                              nrow = nrow(gDiv),
                              ncol = 3))
colnames(dat0D) <- c("species", "geneticDiversity", "type")
for(i in 1:nrow(gDiv)){
  dat0D$species[i] <- gDiv$sp[i]
  dat0D$geneticDiversity[i] <- gDiv$X._0D[i]
  dat0D$type[i] <- host$type[host$sp == gDiv$sp[i]]
}
# plot of 4 fold degenerate sites
plot(x = NULL,
     y = NULL,
     xlim = c(.5,2.5),
     ylim = c(0,max(dat0D$geneticDiversity)),
     xlab = "Host type",
     ylab = "Genetic diversity at zero fold degenerative sites",
     axes = F)
axis(side = 1, at = c(.5,2.5),labels = F)
axis(side = 2)
boxplot(dat0D$geneticDiversity ~ dat0D$type, add = T,axes = F, col = NULL)
points(x =jitter(c(1,2)[dat0D$type + 1]),
       y = dat0D$geneticDiversity,
       cex = 1,
       pch = 16)
axis(side = 1, at = c(1,2), labels = c("Generalist", "Specialist"))
# perform statistical test
t.test <- t.test(dat0D$geneticDiversity ~ dat0D$type)
text(x = 1.5,
     y = 0.0003,
     labels = paste("t test p-value ", round(t.test$p.value, 2), sep = ""),
     pos = 4)
