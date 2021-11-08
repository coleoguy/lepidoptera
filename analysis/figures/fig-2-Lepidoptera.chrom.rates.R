# Terrence Sylvester
# pradakshanas@gmail.com
# August 31, 2021

# load libraries
library(coda)

# load data
amphib.w.xenopus <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv")


# get densities
asc.diff <- density(amphib.w.xenopus$asc1 - amphib.w.xenopus$asc2)
desc.diff <- density(amphib.w.xenopus$desc1 - amphib.w.xenopus$desc2)
pol.diff <- density(amphib.w.xenopus$pol1 - amphib.w.xenopus$pol2)

# get HPD intervals
asc.diff.HPD <- HPDinterval(as.mcmc(amphib.w.xenopus$asc1 - amphib.w.xenopus$asc2))
desc.diff.HPD <- HPDinterval(as.mcmc(amphib.w.xenopus$desc1 - amphib.w.xenopus$desc2))
pol.diff.HPD <- HPDinterval(as.mcmc(amphib.w.xenopus$pol1 - amphib.w.xenopus$pol2))

ylim.max <- max(asc.diff$y,desc.diff$y,pol.diff$y)
ylim.min <-  0 - (max(asc.diff$y,desc.diff$y,pol.diff$y) * 0.1)

# make an empty canvas
plot(x = NULL,
     y = NULL,
     xlim = c(min(asc.diff$x,
                  desc.diff$x,
                  pol.diff$x),
              max(asc.diff$x,
                  desc.diff$x,
                  pol.diff$x)),
     ylim = c(ylim.min,
              max(asc.diff$y,
                  desc.diff$y,
                  pol.diff$y)),
     xlab = expression(paste(Delta, "R"[x])),
     ylab = "Density")
# fissions
polygon(asc.diff, col = rgb(1,0,0,.5),lwd = 2)
segments(x0 = asc.diff.HPD[1,1],
         x1 = asc.diff.HPD[1,2],
         y0 = ylim.min * 2 / 3,
         y1 = ylim.min * 2 / 3,
         col = rgb(1,0,0,.5),
         lwd = 5)
# fusions
polygon(desc.diff, col = rgb(0,1,0,.5), lwd = 2)
segments(x0 = desc.diff.HPD[1,1],
         x1 = desc.diff.HPD[1,2],
         y0 = ylim.min / 3,
         y1 = ylim.min / 3,
         col = rgb(0,1,0,.5),
         lwd = 5)
#poluploidy
polygon(pol.diff, col = rgb(0,0,1,.5), lwd = 2)
segments(x0 = pol.diff.HPD[1,1],
         x1 = pol.diff.HPD[1,2],
         y0 = ylim.min * 3 / 3,
         y1 = ylim.min * 3 / 3,
         col = rgb(0,0,1,.5),
         lwd = 5)
# mark the zero line
abline(v = 0,
       col = "red",
       lty = 2,
       lwd = 2)

# legend
points(x = rep(-0.2, 3),
       y = c(ylim.max, ylim.max - ylim.max * 0.05, ylim.max - ylim.max * 0.1),
       pch = 22,
       bg = c(rgb(1,0,0,.5),
              rgb(0,1,0,.5),
              rgb(0,0,1,.5)),
       col = "black")
# text
text(x = rep(-0.20, 3),
     y = c(ylim.max, ylim.max - ylim.max * 0.05, ylim.max - ylim.max * 0.1), 
     labels = c("Fission", "Fusion", "Polyploidy"),
     pos = 4,
     cex = .9)
# replace chromosome number with branch rates
branch.rates$branch
# isolate tip branches
tip.branches <- matrix(data = NA,
                       nrow = Ntip(tree),
                       ncol = 3)
counter <- 1
for(i in 1:nrow(tree$edge)){
        if(tree$edge[i,2] %in% (1:Ntip(tree))){
                tip.branches[counter,1] <- counter
                tip.branches[counter,2] <- i
                tip.branches[counter,3] <- branch.rates$rate[branch.rates$branch == i]
                counter <- counter + 1
        }
}
MCMC.dat$tipRates <- tip.branches[,3]

# plot tree with chromosome number
plotTree.wBars(tree = tree,
               x = setNames(MCMC.dat$tipRates, MCMC.dat$Species),
               type = "fan",
               col = c("red","blue")[(MCMC.dat$gen.prob + 1)],
               border = NA, 
               lwd = 1,
               scale = 100,
               color = setNames(c("black","black","black",viridis(1,direction = -1)),
                                c(0,1,2,3)),
               method = "plotSimmap")


draw.circle(x = 0,
            y = 0,
            radius = c(100,110,120,130,140),
            nv=100,
            border="gray",
            col=NA,
            lty=2,
            density=NULL,
            angle=45,
            lwd=1)