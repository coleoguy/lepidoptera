# Terrence Sylvester
# 13 October 202
# pradakshanas@g,ail.com

# load libraries
library(phytools)
library(viridis)

# load data
load("../results/stoch.mapping-rates-from-MCMC.RData")

# plot the posterior probability of
DensMap <- densityMap(trees = simmap.MCMC[1:100],
           res = 100,
           ftype = "off",
           lwd = 1.7,
           legend = 50,
           type = "fan", 
           plot = F)

# sort feeding type to match with the order of tip names
feedingType.sorted <- vector(mode = "numeric", length = 441)

for(i in 1:441){
  feedingType.sorted[i] <- host.new$type[host.new$binomial == phy[[1]]$tip.label[i]]
  names(feedingType.sorted)[i] <- phy[[1]]$tip.label[i]
}

plot(DensMap,
     res = 100,
     ftype = "off",
     lwd = 1.7,
     legend = F,
     type = "fan")

# legend
add.color.bar(50,
              DensMap$cols,
              title="Posterior\nporbability (state = 1)",
              lims= c(0,1),
              digits=3,
              prompt=FALSE,
              x=-100, # this is the possition where legend starts. We chose this
              # by using the locator function in R
              y=-91,
              lwd=4,
              fsize=.7,
              subtitle="")

# fix offset
.PlotPhyloEnv$last_plot.phylo$align.tip.label <- T

tiplabels(pch = 16,
          cex = .5, 
          col = c("black", "gray")[feedingType.sorted + 1],
          offset = 1.5)

# plot legend for tip lables
points(x = rep(-100, 2),
       y = c(91,86),
       pch = 16, 
       cex = 1, 
       col = c("black", "gray"))

# add text
text(x = rep(-100, 2),
     y = c(91,86),
     labels = c("Generalist species", "Specialist species"),
     pos = 4,
     cex = .8)

# number of state changes
count <- describe.simmap(simmap.MCMC)$count

# mean transition rates between states
#from generalist to specialist
mean(count[,2])
# from specialist to generalist
mean(count[,3])
