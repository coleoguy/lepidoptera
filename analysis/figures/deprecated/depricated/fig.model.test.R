# Terrence Sylvester
# pradakshanas@gmail.com
# 22 December 2020

# load libraries
library(coda)

# load post burnin from the mcmc and from model test
pb.mcmc <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
pb.model.test <- read.csv("../results/model.test.post.burnin.csv", as.is = T)

# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)

# make a table to get the HPD intervals for each MCMC we did
HPD.tab.fusions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.fissions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.poly <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))

colnames(HPD.tab.fusions) <- colnames(HPD.tab.fissions) <- colnames(HPD.tab.poly) <- c("tree", "HPD-low", "HPD-high", "pass-zero")

# fill in data
for(i in 1:100){
  HPD.tab.fusions$tree[i] <- i
  HPD.tab.fissions$tree[i] <- i
  HPD.tab.poly$tree[i] <- i
  
  HPD.tab.fissions$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$asc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$asc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.fusions$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$desc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$desc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.poly$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$pol1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$pol2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  
  HPD.tab.fissions$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$asc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$asc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.fusions$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$desc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$desc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.poly$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$pol1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$pol2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  
  if(HPD.tab.fusions$`HPD-low`[i] < 0 & HPD.tab.fusions$`HPD-high`[i] > 0){
    HPD.tab.fusions$`pass-zero`[i] <- T  
  }else{
    HPD.tab.fusions$`pass-zero`[i] <- F
  }
  
  if(HPD.tab.fissions$`HPD-low`[i] < 0 & HPD.tab.fissions$`HPD-high`[i] > 0){
    HPD.tab.fissions$`pass-zero`[i] <- T  
  }else{
    HPD.tab.fissions$`pass-zero`[i] <- F
  }
  
  if(HPD.tab.poly$`HPD-low`[i] < 0 & HPD.tab.poly$`HPD-high`[i] > 0){
    HPD.tab.poly$`pass-zero`[i] <- T  
  }else{
    HPD.tab.poly$`pass-zero`[i] <- F
  }
}

# get the rate difference for the mcmc analysis
mcmc.asc.diff <- (pb.mcmc$asc1 - pb.mcmc$asc2)
mcmc.desc.diff <- (pb.mcmc$desc1 - pb.mcmc$desc2)
mcmc.pol.diff <- (pb.mcmc$pol1 - pb.mcmc$pol2)

# get the rate difference for the model testing analysis
mt.asc.diff <- (pb.model.test$asc1 - pb.model.test$asc2)
mt.desc.diff <- (pb.model.test$desc1 - pb.model.test$desc2)
mt.pol.diff <- (pb.model.test$pol1 - pb.model.test$pol2)

# plots from model test
# fission difference
plot(mt.asc.diff,
     cex = .3,
     pch = 16,
     ylab = "Rate difference (MYA)")
segments(x0 = pb.seq,
         x1 = pb.seq,
         y0 = HPD.tab.fissions$`HPD-low`,
         y1 = HPD.tab.fissions$`HPD-high`)
abline(h = 0,
       lwd = 2,
       col = "red")
# fusion difference
plot(mt.desc.diff,
     cex = .3,
     pch = 16,
     ylab = "Rate difference (MYA)")
abline(h = 0,
       lwd = 2,
       col = "red")
# polyploidy difference
plot(mt.pol.diff,
     cex = .3,
     pch = 16,
     ylab = "Rate difference (MYA)")
abline(h = 0,
       lwd = 2,
       col = "red")
segments(x0 = c(1,100,200),
         x1 = c(1,100,200),
         y0 = 0, y = .04)



par(mfcol = c(1,3))
# make an empty canvas for fissions
plot(x = NULL,
     y = NULL,
     xlim = c(-.4, .4),
     ylim = c(0,max(density(mt.asc.diff)$y, density(mcmc.asc.diff)$y)),
     xlab = "Difference of the rates of fission (MYA)",
     ylab = "Density",
     cex.lab = 1.5,
     cex.axis = 1.5)
  mtext(text = "A",
        line = 0,
        outer = F,
        side = 3,
        adj = 0,
        cex = 1.5)
polygon(density(mt.asc.diff),lwd = 2, col = rgb(1,0,0,.5))
polygon(density(mcmc.asc.diff), lwd = 2, col = rgb(0,0,1,.5))

abline(v = mean(mcmc.asc.diff), col = "red", lwd = 3)
# make an empty canvas for fusions
plot(x = NULL,
     y = NULL,
     xlim = c(-.7, .7),
     ylim = c(0,max(mt.desc.diff$y)),
     xlab = "Difference of the rates of fusion (MYA)",
     ylab = "Density",
     cex.lab = 1.5,
     cex.axis = 1.5)
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1.5)
polygon(mt.desc.diff,lwd = 2)
abline(v = mean(mcmc.desc.diff), col = "red", lwd = 3)
# make an empty canvas for polyploidy
plot(x = NULL,
     y = NULL,
     xlim = c(-.08, .08),
     ylim = c(-0,max(mt.pol.diff$y)),
     xlab = "Difference of the rates of polyploidy (MYA)",
     ylab = "Density",
     cex.lab = 1.5,
     cex.axis = 1.5)
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1.5)
polygon(mt.pol.diff,lwd = 2)
abline(v = mean(mcmc.pol.diff), col = "red", lwd = 3)



