library(phytools)

tree <- read.newick("../Hardy-Otto-Data/Papilionoidea-BiSSE/papilionoidea-10.trees")
dat <- read.csv("../Hardy-Otto-Data/Papilionoidea-BiSSE/papilionoidea-bisse-data.txt")

tree <- tree[[1]]

plot(tree,
     show.tip.label = F,
     type = "f")


tiplabels(pch = 16,
          col = c("red", "black")[dat$hosts + 1],
          offset = 0.03,
          cex = .5)

legend(x = 100.3906,
       y = 94.24288,
       legend = c("Specialist", "Generalist"),
       pch = 16,
       col = c("black", "red"),
       bty = "n")

# data oerlap

keep.tips <- unique(chroms$species)

new.tree <-  drop.tip(phy = tree, 
                      tip = keep.tips)

new.tree <- drop.tip(phy = tree,
                     tip = new.tree$tip.label)

new.dat <- dat[dat$sp %in% new.tree$tip.label,]

plot(new.tree, show.tip.label = F,
     type = "f", edge.width = 2)

tiplabels(pch = 16,
          col = c("red", "black")[new.dat$hosts + 1],
          offset = 0.05,
          cex = 1)

legend(x = 94,
       y = 94,
       legend = c("Specialist", "Generalist"),
       pch = 16,
       col = c("black", "red"),
       bty = "n")
