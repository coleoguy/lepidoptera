# load libraries
library(phytools)
library(viridis)

# load data
load("../results/9b.branch.specific.rates.Lepidoptera.with.feeding.RData")

# process branch rates
# Assign a rate class for each branch depending on the branch specific rates
for(i in 1:nrow(branch.rates)){
        if(branch.rates$rate[i] <= quantile(x = branch.rates$rate, 0.25)){
                branch.rates$rateClass[i] <- 1
        }
        if(branch.rates$rate[i] > quantile(x = branch.rates$rate, 0.25) & branch.rates$rate[i] < quantile(x = branch.rates$rate, 0.9)){
                branch.rates$rateClass[i] <- 2
        }
        if(branch.rates$rate[i] >= quantile(x = branch.rates$rate, 0.9)){
                branch.rates$rateClass[i] <- 3
        }
}
# colour the tip branches only
branch.rates$rateClass[!(tree$edge[,2] %in% (1:441))] <- 0

# make a maps object and add it to tree in order to colour branches
maps <- vector(mode = "list", length = Nedge(tree))
for(i in 1:length(maps)){
        maps[[i]] <- tree$edge.length[i]
        names(maps[[i]]) <- branch.rates$rateClass[i]
}
tree$maps <- maps
# plot tree with chromosome number
plotTree.wBars(tree = tree,
               x = setNames(MCMC.dat$Chroms, MCMC.dat$Species),
               type = "fan",
               col = c("red","blue")[(MCMC.dat$gen.prob + 1)],
               border = NA, 
               lwd = 1,
               color = setNames(c("black",viridis(3)),
                                c(0,1,2,3)),
               method = "plotSimmap")

# # plot tiplables
# # get the overlap between the species and chromosome number data sets
# obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
# n<-Ntip(tree)
# col.bin<-setNames(c("red","blue"),levels(as.factor(MCMC.dat$gen.prob)))
# 
# offset <- 0.03
# par(lend=3)
# for(i in 1:Ntip(tree)){
#         if(max(branching.times(tree)) < 1){
#                 cc<-if(obj$xx[i]>0) 5 else -5
#                 th<-atan(obj$yy[i]/obj$xx[i])
#                 segments(obj$xx[i] + (obj$xx[i] * offset),
#                          obj$yy[i] + (obj$yy[i] * offset),
#                          obj$xx[i]+(cc/tree.depth)*cos(th) + (obj$xx[i] * offset),
#                          obj$yy[i]+(cc/tree.depth)*sin(th) + (obj$yy[i] * offset),
#                          lwd=2,
#                          col=col.bin[MCMC.dat$gen.prob[MCMC.dat$Species == tree$tip.label[i]] + 1])
#                 
#         }
#         else{
#                 cc<-if(obj$xx[i]>0) 5 else -5
#                 th<-atan(obj$yy[i]/obj$xx[i])
#                 segments(obj$xx[i] + (obj$xx[i] * offset),
#                          obj$yy[i] + (obj$yy[i] * offset),
#                          obj$xx[i]+cc*cos(th) + (obj$xx[i] * offset),
#                          obj$yy[i]+cc*sin(th) + (obj$yy[i] * offset),
#                          lwd=2,
#                          col=col.bin[MCMC.dat$gen.prob[MCMC.dat$Species == tree$tip.label[i]] + 1])
#         }
# }
# plot the legend
text(x = 0.75,
     y = 1.05 + 0.25,
     labels = "Rate class",
     pos = 4,
     cex = 0.8)
points(x = rep(.8,3),
       y = c(1.0,0.95,0.90)+ 0.25,
       pch = 16,
       col = viridis(3))
text(x = rep(.8,3),
     y = c(1.0,0.95,0.90)+ 0.25,
     labels = c(paste("Low (<= ",
                      round(quantile(x = branch.rates$rate, 0.25)/ tree.depth,2),
                      ")"),
                paste("Medium (< ",
                      round(quantile(x = branch.rates$rate, 0.25)/ tree.depth,2),
                      "and",
                      round(quantile(x = branch.rates$rate, 0.9)/ tree.depth,2),
                      ">)"),
                paste("High (>= ",
                      round(quantile(x = branch.rates$rate, 0.9)/ tree.depth,2)
                      ,")")),
     pos = 4, cex = .7)
points(x = rep(.8,2),
       y = c(0.79,0.74)+ 0.23,
       pch = 16,
       col = c("blue","red"))
text(x = 0.75,
     y = 0.84+ 0.23,
     labels = "Binary trait",
     pos = 4, cex = .8)
text(x = rep(.8,2),
     y = c(0.79,0.74) + 0.23,
     labels = c("Generalists", "Specialists"),
     pos = 4,
     cex = 0.7)



