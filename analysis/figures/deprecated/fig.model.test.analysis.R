# Terrence Sylvester
# pradakshanas@gmail.com
# 22 December 2020

# load libraries
library(coda)

# load post burnin from the mcmc and from model test
pb.model.test <- read.csv("../results/model.test.post.burnin.csv", as.is = T)
pb.mcmc <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)

# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)

# make a table to get the HPD intervals for each MCMC we did
HPD.tab.fusions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.fissions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.poly <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))

colnames(HPD.tab.fusions) <- colnames(HPD.tab.fissions) <- colnames(HPD.tab.poly) <- c("tree", "HPD-low", "HPD-high", "pass-zero")

# make a table to get the distribution of mean delta r statistic
Delta.R <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
colnames(Delta.R) <- c("Sim", "asc", "desc", "pol")



# fill in data
for(i in 1:100){
  HPD.tab.fusions$tree[i] <- i
  HPD.tab.fissions$tree[i] <- i
  HPD.tab.poly$tree[i] <- i
  
  HPD.tab.fusions$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$asc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$asc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.fissions$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$desc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$desc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.poly$`HPD-low`[i] <- HPDinterval(as.mcmc(pb.model.test$pol1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$pol2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  
  HPD.tab.fusions$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$asc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$asc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.fissions$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$desc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$desc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.poly$`HPD-high`[i] <- HPDinterval(as.mcmc(pb.model.test$pol1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$pol2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  
  Delta.R$Sim[i] <- i
  Delta.R$asc[i] <- mean(abs(pb.model.test$asc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$asc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$desc[i] <- mean(abs(pb.model.test$desc1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$desc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$pol[i] <- mean(abs(pb.model.test$pol1[pb.seq[i]:(pb.seq[i]+49)] - pb.model.test$pol2[pb.seq[i]:(pb.seq[i]+49)]))
  
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

sum(HPD.tab.fusions$`pass-zero`)
sum(HPD.tab.fusions$`pass-zero`)
sum(HPD.tab.poly$`pass-zero`)

par(mfcol = c(1,3))

# fission
plot(density(Delta.R$asc), 
     xlab = expression(paste("| ",Delta, "R"[fission]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(pb.mcmc$asc1 - pb.mcmc$asc2),
       col = "red",
       lwd = 2)
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# fusion
plot(density(Delta.R$desc), 
     xlab = expression(paste("| ",Delta, "R"[fusion]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(pb.mcmc$desc1 - pb.mcmc$desc2),
       col = "red",
       lwd = 2)
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# polyploidy
plot(density(Delta.R$pol), 
     xlab = expression(paste("| ",Delta, "R"[polyploidy]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(pb.mcmc$pol1 - pb.mcmc$pol2),
       col = "red",
       lwd = 2)
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)

sum(Delta.R$asc > mean(pb.mcmc$asc1 - pb.mcmc$asc2))
sum(Delta.R$desc > mean(pb.mcmc$desc1 - pb.mcmc$desc2))
sum(Delta.R$pol > mean(pb.mcmc$pol1 - pb.mcmc$pol2))
