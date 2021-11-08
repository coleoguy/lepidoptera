# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# stochastic mapping

# load libraries
library(phytools)

# read in trees
trees <- read.nexus("../data/trees/heliconiini-100.nex")

# read in hosts data
hosts <- read.csv("../data/hosts/heliconiini-bisse.csv", as.is = T)

simmaps <- simmap.summaries <- vector("list", length = 100)
names(simmaps) <- names(simmap.summaries) <- paste("tree", 1:100, sep = "")

for(i in 1:100){
  nsim <- 100
  
  # sample a single tree
  tree <- trees[[i]]
  
  # keep tips for which we have data
  tree <- keep.tip(phy = tree, tip = hosts$sp)
  
  # scale tree to unit length
  tree.depth <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth
  
  # track progress
  print(paste("working on tree ", i))

  # estimate ancestral states
  chars <- hosts$range
  names(chars) <- hosts$sp
  
  # perform stochastic mapping
  simmaps[[i]] <- make.simmap(tree = tree,
                              x = chars,
                              nsim = nsim,
                              model = "ARD")
  
  # get summaries
  simmap.summaries[[i]] <- describe.simmap(simmaps[[i]])
  
}

