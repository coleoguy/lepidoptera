# Terrence Sylvester
# 13 October 202
# pradakshanas@g,ail.com

# load libraries
library(phytools)
library(viridis)

# load data
load("../results/stoch.mapping.RData")

# plot the posterior probability of
DensMap <- densityMap(trees = simmap$tree.1,
           res = 100,
           ftype = "off",
           lwd = 1.7,
           legend = 50,
           type = "fan", 
           plot = F)

# sort feeding type to match with the order of tip names
feedingType.sorted <- vector(mode = "numeric", length = 1039)

for(i in 1:1039){
  feedingType.sorted[i] <- hosts$type[hosts$spNames == trees[[1]]$tip.label[i]]
  names(feedingType.sorted)[i] <- trees[[1]]$tip.label[i]
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
          cex = .3, 
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
count <- matrix(data = NA,
                nrow = 0,
                ncol = 3)

for(i in 1:10){
count <- rbind(describe.simmap(simmap[[i]])$count, count)
}

# mean transition rates between states
#from generalist to specialist
mean(count[,2])
# from specialist to generalist
mean(count[,3])
