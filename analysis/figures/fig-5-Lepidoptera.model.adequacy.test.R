# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# read in data
load("../results/6a.model.adequacy.test.Lepidoptera.rate.analysis.feeding.RData")

# get the proportion of non direct developing taxa 
par(mfcol = c(1,3))
prop.Gen <- c()
for(i in 1:100){
  for(j in 1:100){
    prop.Gen <-  c(prop.Gen, sum(results.final[[i]]$states[[j]])/length(results.final[[i]]$states[[j]]))
  }
}
#plot the proportion of non direct development
plot(density(prop.Gen),
     main = "",
     xlab = "Proportion of generalists",
     lwd = 2)
abline(v = sum(MCMC.dat$gen.prob)/nrow(MCMC.dat),
       col = "red",
       lwd = 2)
#plot title
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# observed coefficient of variance in chromosome number
obs.coef.var <- sd(MCMC.dat$Chroms) / mean(MCMC.dat$Chroms)
# simulated coefficient of variance in chromosome number
sim.coef.var <- c()
for(i in 1:100){
  for(j in 1:100){
    sim.coef.var <- c(sim.coef.var,(sd(results.final[[i]]$chroms[[j]])/mean(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.coef.var),
     xlab = "Coefficient of variance in chromosome number",
     main = "",
     lwd = 2)
abline(v = obs.coef.var,
       col = "red",
       lwd = 2)
#plot title
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# observed variance in chromosome number
obs.var <- var(MCMC.dat$Chroms)
# simulated variance in chromosome number
sim.var <- c()
for(i in 1:100){
  for(j in 1:100){
    sim.var <- c(sim.var,(var(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.var),
     xlab = "Variance in chromosome number",
     main = "",
     lwd = 2)
abline(v = obs.var,
       col = "red",
       lwd = 2)
#plot title
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)

