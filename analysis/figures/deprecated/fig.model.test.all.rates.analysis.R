# load libraries
library(coda)

# load post burnin from the mcmc and from model test
pb.simulated.test <- read.csv("../results/model.test.post.burnin.csv", as.is = T)
pb.mcmc <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
# get the sum of all rates at each state
emperical.state.01.rates <- pb.mcmc$asc1 + pb.mcmc$desc1 + pb.mcmc$pol1
emperical.state.02.rates <- pb.mcmc$asc2 + pb.mcmc$desc2 + pb.mcmc$pol2
simulated.state.01.rates <- pb.simulated.test$asc1 + pb.simulated.test$desc1 + pb.simulated.test$pol1
simulated.state.02.rates <- pb.simulated.test$asc2 + pb.simulated.test$desc2 + pb.simulated.test$pol2
# compare these rates
emperical.rate.diff <- emperical.state.01.rates - emperical.state.02.rates
simulated.rate.diff <- abs(simulated.state.01.rates - simulated.state.02.rates)
# density for plotting
emperical.rate.density <- density(emperical.rate.diff)
simulated.rate.density <- density(simulated.rate.diff)
# HPD interval
emperical.rate.HPD <- HPDinterval(as.mcmc(emperical.rate.diff))
simulated.rate.HPD <- HPDinterval(as.mcmc(simulated.rate.diff))
# plot
plot(x = NULL,
     y = NULL,
     xlim = c(0,1),
     ylim = c(-1,10),
     xlab = expression(paste(Delta, "R"[All-eperical], " vs. ", Delta, "R"[All-simulated], " (MYA)")),
     main = "")
polygon(emperical.rate.density,
        col = rgb(0,0,1,.5))
polygon(simulated.rate.density,
        col = rgb(1,0,0,.5))

segments(x0 = emperical.rate.HPD[1],
         x1 = emperical.rate.HPD[2],
         y0 = -.5,
         y1 = -.5,
         lwd = 3,
         col = "Blue")
segments(x0 = simulated.rate.HPD[1],
         x1 = simulated.rate.HPD[2],
         y0 = -1,
         y1 = -1,
         lwd = 3,
         col = "Red")
abline(v = 0,
       col = "red",
       lty = 2)
# legend
points(x = rep(0.8, 2),
       y = c(10,9),
       pch = 22,
       bg = c(rgb(0,0,1,.5),
              rgb(1,0,0,.5)),
       col = "black",
       cex = 2)
# text
text(x = rep(0.8, 2),
     y = c(10,9), 
     labels = c("Emperical dataset", "Simulated dataset"),
     pos = 4,
     cex = .9)
