# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# stochastic mapping

# load libraries
library(diversitree)
library(phytools)

# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

# read in hosts data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)

simmaps <- simmap.summaries <- vector("list", length = 10)

for(i in 1:10){
  iter <- 1000
  
  # rename tree tips to match with trait data set
  tree <- trees[[i]]
  # when visually axamining the tree, these trees look ultrametric.(they shoud be) 
  # however the when making the likelihood functions if does not recognize these
  # trees as being ultrametric. therefore we will force the trees to be
  # ultrametric
  tree <- force.ultrametric(tree = tree)
  
  # scale tree to unit length
  tree.depth <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth
  
  # track progress
  print(paste("working on tree ", i))

  # estimate ancestral states
  chars <- hosts$hosts
  names(chars) <- hosts$sp
  
  simmaps[[i]] <- make.simmap(tree = tree,
                              x = chars,
                              nsim = 100,
                              model = "ARD")
  
  simmap.summaries[[i]] <- describe.simmap(simmaps[[i]])
  
}

