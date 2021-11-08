# Terrence Sylvester
# pradakshanas@gmail.com
# May 07, 2020
# papilionoidea-mcmc-mkn

# load libraries
library(coda)

# load data
dat <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
# get the rate difference between the two states
# state 1 is specialists and state 2 is generalist
# we will get the difference by looking at state 1 - state 2
asc.diff <- density(dat$asc1 - dat$asc2)
desc.diff <- density(dat$desc1 - dat$desc2)
pol.diff <- density(dat$pol1 - dat$pol2)
# get the HPD values
asc.diff.HPD <- HPDinterval(as.mcmc(dat$asc1 - dat$asc2))
desc.diff.HPD <- HPDinterval(as.mcmc(dat$desc1 - dat$desc2))
pol.diff.HPD <- HPDinterval(as.mcmc(dat$pol1 - dat$pol2))

# make an empty canvas
plot(x = NULL,
     y = NULL,
     xlim = c(0, .4),
     ylim = c(-6,70),
     xlab = "Difference of the rates of \nchromosome number evolution (MYA)",
     ylab = "Density")
# fissions
polygon(asc.diff, col = rgb(1,0,0,.5),lwd = 2)
segments(x0 = asc.diff.HPD[1,1],
         x1 = asc.diff.HPD[1,2],
         y0 = -4,
         y1 = -4,
         col = rgb(1,0,0,.5),
         lwd = 5)
# fusions
polygon(desc.diff, col = rgb(0,1,0,.5), lwd = 2)
segments(x0 = desc.diff.HPD[1,1],
         x1 = desc.diff.HPD[1,2],
         y0 = -6,
         y1 = -6,
         col = rgb(0,1,0,.5),
         lwd = 5)
#poluploidy
polygon(pol.diff, col = rgb(0,0,1,.5), lwd = 2)
segments(x0 = pol.diff.HPD[1,1],
         x1 = pol.diff.HPD[1,2],
         y0 = -2,
         y1 = -2,
         col = rgb(0,0,1,.5),
         lwd = 5)
# mark the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)
# legend
points(x = rep(0.35, 3),
       y = c(70,67,64),
       pch = 22,
       bg = c(rgb(1,0,0,.5),
               rgb(0,1,0,.5),
               rgb(0,0,1,.5)),
       col = "black")
# text
text(x = rep(0.35, 3),
     y = c(70,67,64), 
     labels = c("Fission", "Fusion", "Polyploidy"),
     pos = 4,
     cex = .9)
