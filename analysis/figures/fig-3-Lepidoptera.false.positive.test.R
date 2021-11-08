# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# load libraries
library(coda)

# load data
emp.dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv")
sim.dat <- read.csv("../results/4a.false.positive.test.Lepidoptera.rate.analysis.feeding.csv")

# get densities
emp.asc.diff <- density(emp.dat$asc1 - emp.dat$asc2)
emp.desc.diff <- density(emp.dat$desc1 - emp.dat$desc2)
emp.pol.diff <- density(emp.dat$pol1 - emp.dat$pol2)

sim.asc.diff <- density(sim.dat$asc1 - sim.dat$asc2)
sim.desc.diff <- density(sim.dat$desc1 - sim.dat$desc2)
sim.pol.diff <- density(sim.dat$pol1 - sim.dat$pol2)


# get HPD intervals
emp.asc.diff.HPD <- HPDinterval(as.mcmc(emp.dat$asc1 - emp.dat$asc2))
emp.desc.diff.HPD <- HPDinterval(as.mcmc(emp.dat$desc1 - emp.dat$desc2))
emp.pol.diff.HPD <- HPDinterval(as.mcmc(emp.dat$pol1 - emp.dat$pol2))

sim.asc.diff.HPD <- HPDinterval(as.mcmc(sim.dat$asc1 - sim.dat$asc2))
sim.desc.diff.HPD <- HPDinterval(as.mcmc(sim.dat$desc1 - sim.dat$desc2))
sim.pol.diff.HPD <- HPDinterval(as.mcmc(sim.dat$pol1 - sim.dat$pol2))


par(mfcol = c(1,3))

# compare fission rates
ylim.max <- max(emp.asc.diff$y,sim.asc.diff$y)
ylim.min <-  0 - (max(emp.asc.diff$y,sim.asc.diff$y) * 0.1)

# make an empty canvas
plot(x = NULL,
     y = NULL,
     xlim = c(min(emp.asc.diff$x,
                  sim.asc.diff$x),
              max(emp.asc.diff$x,
                  sim.asc.diff$x)),
     ylim = c(ylim.min,
              ylim.max),
     xlab = expression(paste(Delta, "R"[X])),
     ylab = "Density")
#plot title
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# Emperical fission
polygon(emp.asc.diff, col = rgb(1,0,0,.5),lwd = 2)
segments(x0 = emp.asc.diff.HPD[1,1],
         x1 = emp.asc.diff.HPD[1,2],
         y0 = ylim.min * 2 / 3,
         y1 = ylim.min * 2 / 3,
         col = rgb(1,0,0,.5),
         lwd = 2)
# simulated fission
polygon(sim.asc.diff, col = rgb(.1,.1,.1,.5), lwd = 2)
segments(x0 = sim.asc.diff.HPD[1,1],
         x1 = sim.asc.diff.HPD[1,2],
         y0 = ylim.min / 3,
         y1 = ylim.min / 3,
         col = rgb(.1,.1,.1,.5),
         lwd = 2)
# mark the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# compare fusion rates
ylim.max <- max(emp.desc.diff$y,sim.desc.diff$y)
ylim.min <-  0 - (max(emp.desc.diff$y,sim.desc.diff$y) * 0.1)

# make an empty canvas
plot(x = NULL,
     y = NULL,
     xlim = c(min(emp.desc.diff$x,
                  sim.desc.diff$x),
              max(emp.desc.diff$x,
                  sim.desc.diff$x)),
     ylim = c(ylim.min,
              ylim.max),
     xlab = expression(paste(Delta, "R"[X])),
     ylab = "Density")
#plot title
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# Emperical fission
polygon(emp.desc.diff, col = rgb(0,1,0,.5),lwd = 2)
segments(x0 = emp.desc.diff.HPD[1,1],
         x1 = emp.desc.diff.HPD[1,2],
         y0 = ylim.min * 2 / 3,
         y1 = ylim.min * 2 / 3,
         col = rgb(0,1,0,.5),
         lwd = 2)
# simulated fission
polygon(sim.desc.diff, col = rgb(.1,.1,.1,.5), lwd = 2)
segments(x0 = sim.desc.diff.HPD[1,1],
         x1 = sim.desc.diff.HPD[1,2],
         y0 = ylim.min / 3,
         y1 = ylim.min / 3,
         col = rgb(.1,.1,.1,.5),
         lwd = 2)
# mark the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# compare polyploidy  rates
ylim.max <- max(emp.pol.diff$y,sim.pol.diff$y)
ylim.min <-  0 - (max(emp.pol.diff$y,sim.pol.diff$y) * 0.1)

# make an empty canvas
plot(x = NULL,
     y = NULL,
     xlim = c(min(emp.pol.diff$x,
                  sim.pol.diff$x),
              max(emp.pol.diff$x,
                  sim.pol.diff$x)),
     ylim = c(ylim.min,
              ylim.max),
     xlab = expression(paste(Delta, "R"[X])),
     ylab = "Density")
#plot title
mtext(text = "C",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# Emperical fission
polygon(emp.pol.diff, col = rgb(0,0,1,.5),lwd = 2)
segments(x0 = emp.pol.diff.HPD[1,1],
         x1 = emp.pol.diff.HPD[1,2],
         y0 = ylim.min * 2 / 3,
         y1 = ylim.min * 2 / 3,
         col = rgb(0,0,1,.5),
         lwd = 2)
# simulated fission
polygon(sim.pol.diff, col =rgb(.1,.1,.1,.5), lwd = 2)
segments(x0 = sim.pol.diff.HPD[1,1],
         x1 = sim.pol.diff.HPD[1,2],
         y0 = ylim.min / 3,
         y1 = ylim.min / 3,
         col = rgb(.1,.1,.1,.5),
         lwd = 2)
# mark the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# legend
points(x = rep(0.01,4),
       y = c(ylim.max, 
             ylim.max - ylim.max * 0.05,
             ylim.max - ylim.max * 0.1,
             ylim.max - ylim.max * 0.15),
       pch = 22,
       bg = c(rgb(1,0,0,.5),
              rgb(0,1,0,.5),
              rgb(0,0,1,.5),
              rgb(.1,.1,.1,.5)),
       col = "black")
# text
text(x = rep(0.01,4),
     y = c(ylim.max, 
           ylim.max - ylim.max * 0.05,
           ylim.max - ylim.max * 0.1,
           ylim.max - ylim.max * 0.15), 
     labels = c("Fusion",
                "Fission",
                "Polyploidy",
                "Simulated dataset"),
     pos = 4,
     cex = .9)

