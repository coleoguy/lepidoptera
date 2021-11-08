# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# load libraries
library(coda)

# load data
emp.dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv")
sim.dat <- read.csv("../results/4a.false.positive.test.Lepidoptera.rate.analysis.feeding.proc.csv")

# make a table to get the HPD intervals for each MCMC we did
HPD.tab.fusions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.fissions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.poly <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
# get column names
colnames(HPD.tab.fusions) <- colnames(HPD.tab.fissions) <- colnames(HPD.tab.poly) <- c("tree", "HPD-low", "HPD-high", "pass-zero")
# make a table to get the distribution of mean delta r statistic
Delta.R <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
colnames(Delta.R) <- c("Sim", "asc", "desc", "pol")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# fill in data
for(i in 1:100){
  # get the working tree
  HPD.tab.fusions$tree[i] <- i
  HPD.tab.fissions$tree[i] <- i
  HPD.tab.poly$tree[i] <- i
  # get the lower limit of HPD
  HPD.tab.fusions$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.fissions$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.poly$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  # get the higher limit of HPD
  HPD.tab.fusions$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.fissions$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.poly$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  # get the mean of the delta R statistic for each run
  Delta.R$Sim[i] <- i
  Delta.R$asc[i] <- mean(abs(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$desc[i] <- mean(abs(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$pol[i] <- mean(abs(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))
  
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
sum(HPD.tab.fissions$`pass-zero`)
sum(HPD.tab.poly$`pass-zero`)

# plot the results
par(mfcol = c(1,3))
# fission
plot(density(Delta.R$asc), 
     xlab = expression(paste("| ",Delta, "R"[fission]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$asc1 - emp.dat$asc2)),
       col = "red",
       lwd = 2)
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
text(x = min(density(Delta.R$asc)$x),
     y = max(density(Delta.R$asc)$y),
     label = paste("p-value",sum(Delta.R$asc > mean(abs(emp.dat$asc1 - emp.dat$asc2))) / 100),
     pos = 4,
     cex = 1.2) 
# fusion
plot(density(Delta.R$desc), 
     xlab = expression(paste("| ",Delta, "R"[fusion]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$desc1 - emp.dat$desc2)),
       col = "red",
       lwd = 2)
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
text(x = max(density(Delta.R$desc)$x),
     y = max(density(Delta.R$desc)$y),
     label = paste("p-value",sum(Delta.R$desc > mean(abs(emp.dat$desc1 - emp.dat$desc2))) / 100),
     pos = 2,
     cex = 1.2) 

# polyploidy
plot(density(Delta.R$pol), 
     xlab = expression(paste("| ",Delta, "R"[polyploidy]," |")),
     main = "",
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$pol1 - emp.dat$pol2)),
       col = "red",
       lwd = 2)
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
text(x = min(density(Delta.R$pol)$x),
     y = max(density(Delta.R$pol)$y),
     label = paste("p-value",sum(Delta.R$pol > mean(abs(emp.dat$pol1 - emp.dat$pol2))) / 100),
     pos = 4,
     cex = 1.2) 


# histograms
# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# load libraries
library(coda)

# load data
emp.dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv")
sim.dat <- read.csv("../results/4a.false.positive.test.Lepidoptera.rate.analysis.feeding.proc.csv")

# make a table to get the HPD intervals for each MCMC we did
HPD.tab.fusions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.fissions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.poly <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
# get column names
colnames(HPD.tab.fusions) <- colnames(HPD.tab.fissions) <- colnames(HPD.tab.poly) <- c("tree", "HPD-low", "HPD-high", "pass-zero")
# make a table to get the distribution of mean delta r statistic
Delta.R <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
colnames(Delta.R) <- c("Sim", "asc", "desc", "pol")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# fill in data
for(i in 1:100){
  # get the working tree
  HPD.tab.fusions$tree[i] <- i
  HPD.tab.fissions$tree[i] <- i
  HPD.tab.poly$tree[i] <- i
  # get the lower limit of HPD
  HPD.tab.fusions$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.fissions$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.poly$`HPD-low`[i] <- HPDinterval(as.mcmc(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  # get the higher limit of HPD
  HPD.tab.fusions$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.fissions$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.poly$`HPD-high`[i] <- HPDinterval(as.mcmc(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  # get the mean of the delta R statistic for each run
  Delta.R$Sim[i] <- i
  Delta.R$asc[i] <- mean(abs(sim.dat$asc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$asc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$desc[i] <- mean(abs(sim.dat$desc1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$desc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$pol[i] <- mean(abs(sim.dat$pol1[pb.seq[i]:(pb.seq[i]+49)] - sim.dat$pol2[pb.seq[i]:(pb.seq[i]+49)]))
  
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
sum(HPD.tab.fissions$`pass-zero`)
sum(HPD.tab.poly$`pass-zero`)

# plot the results
par(mfcol = c(1,3))
# fission
hist(Delta.R$asc, 
     xlab = expression(paste("| ",Delta, "R"[fission]," |")),
     main = "",
     breaks = 25,
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$asc1 - emp.dat$asc2)),
       col = "red",
       lwd = 2)
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# text(x = min(density(Delta.R$asc)$x),
#      y = max(density(Delta.R$asc)$y),
#      label = paste("p-value",sum(Delta.R$asc > mean(abs(emp.dat$asc1 - emp.dat$asc2))) / 100),
#      pos = 4,
#      cex = 1.2) 
# fusion
hist(Delta.R$desc, 
     xlab = expression(paste("| ",Delta, "R"[fusion]," |")),
     main = "",
     breaks = 25,
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$desc1 - emp.dat$desc2)),
       col = "red",
       lwd = 2)
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# text(x = max(density(Delta.R$desc)$x),
#      y = max(density(Delta.R$desc)$y),
#      label = paste("p-value",sum(Delta.R$desc > mean(abs(emp.dat$desc1 - emp.dat$desc2))) / 100),
#      pos = 2,
#      cex = 1.2) 

# polyploidy
hist(Delta.R$pol, 
     xlab = expression(paste("| ",Delta, "R"[polyploidy]," |")),
     main = "",
     breaks = 25,
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(emp.dat$pol1 - emp.dat$pol2)),
       col = "red",
       lwd = 2)
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# text(x = min(density(Delta.R$pol)$x),
#      y = max(density(Delta.R$pol)$y),
#      label = paste("p-value",sum(Delta.R$pol > mean(abs(emp.dat$pol1 - emp.dat$pol2))) / 100),
#      pos = 4,
#      cex = 1.2) 



