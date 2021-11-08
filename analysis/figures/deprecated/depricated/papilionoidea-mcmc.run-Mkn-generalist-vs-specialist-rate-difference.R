# Terrence Sylvester
# pradakshanas@gmail.com
# May 07, 2020
# papilionoidea-mcmc-mkn

# load libraries
library(coda)

# load data
dat <- read.csv("../results/papilionoidea-mcmc.run-Mkn-post-burn-in.csv",
                as.is = T)

# get the densities fro each parameter
# generalist.asc <- density(dat$asc1)
# generalist.desc <- density(dat$desc1)
# generalist.pol <- density(dat$pol1)
# specialist.asc <- density(dat$asc2)
# specialist.desc <- density(dat$desc2)
# specialist.pol <- density(dat$pol2)
# generalist.to.specialist <- density(dat$tran12)
# specialist.to.generalist <- density(dat$tran21)

# get the rate difference between the two states
# state 1 is generalists and state 2 is specialists
# we will get the difference by looking at state 1 - state 2

asc.diff <- density(dat$asc1 - dat$asc2)
desc.diff <- density(dat$desc1 - dat$desc2)
pol.diff <- density(dat$pol1 - dat$pol2)

# get the HPD values
asc.diff.HPD <- HPDinterval(as.mcmc(dat$asc1 - dat$asc2))
desc.diff.HPD <- HPDinterval(as.mcmc(dat$desc1 - dat$desc2))
pol.diff.HPD <- HPDinterval(as.mcmc(dat$pol1 - dat$pol2))

# devide plot window
par(mfcol = c(1,3))

##### chromosome fission in generalists vs specialists ####
plot(x = NULL,
     y = NULL,
     xlim = c(-0.22, 0.1),
     ylim = c(-.8,8),
     xlab = "Difference in rates (MYA)",
     ylab = "Probability density",
     main = "",
     cex.lab = 1.5,
     cex.axis = 1.2)

polygon(asc.diff,
        col = rgb(0,0,1,.4),
        border = rgb(0,0,1,1))

# plot the HPD interval
segments(x0 = asc.diff.HPD[1,1],
         x1 = asc.diff.HPD[1,2],
         y0 = -0.4,
         y1 = -0.4,
         col = rgb(0,0,1,1),
         lwd = 5)

# plot the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# plot lable
mtext(text  = "Fission",
      at = c(-0.2,37.8),
      cex = 2)

##### chromosome fusion in generalists vs specialists ####
plot(x = NULL,
     y = NULL,
     xlim = c(-0.23, 0.05),
     ylim = c(-1.5,15),
     xlab = "Difference in rates (MYA)",
     ylab = "Probability density",
     main = "",
     cex.lab = 1.5,
     cex.axis = 1.2)

polygon(desc.diff,
        col = rgb(0,0,1,.4),
        border = rgb(0,0,1,1))

# plot the HPD interval
segments(x0 = desc.diff.HPD[1,1],
         x1 = desc.diff.HPD[1,2],
         y0 = -.75,
         y1 = -.75,
         col = rgb(0,0,1,1),
         lwd = 5)

# plot the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# plot lable
mtext(text  = "Fusion",
      at = c(-0.2,37.8),
      cex = 2)

#### chromosome fission in generalists vs specialists ####
plot(x = NULL,
     y = NULL,
     xlim = c(-0.09, 0.02),
     ylim = c(-3,30),
     xlab = "Difference in rates (MYA)",
     ylab = "Probability density",
     main = "",
     cex.lab = 1.5,
     cex.axis = 1.2)

polygon(pol.diff,
        col = rgb(0,0,1,.4),
        border = rgb(0,0,1,1))

# plot the HPD interval
segments(x0 = pol.diff.HPD[1,1],
         x1 = pol.diff.HPD[1,2],
         y0 = -1.5,
         y1 = -1.5,
         col = rgb(0,0,1,1),
         lwd = 5)

# plot the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# plot lable
mtext(text  = "Polyploidy",
      at = c(-0.08,97.2),
      cex = 2)

#####

# add title and other notations
# title(main = "Differences in the rates of chromosome number evolution between generalists and specialist species", 
#       outer = T,
#       line = -2,
#       cex.main = 1.5)
# 


