# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# stochastic mapping

# load libraries
library(diversitree)
library(phytools)

# load helper functions
source("helper.functions.R")

# read in chromosome data
dat.initial <- read.csv("../data/chroms/lep.chroms.csv", as.is = T)

# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

# read in hosts data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)

simmaps <- simmap.summaries <- vector("list", length = 100)

for(i in 1:10){
  iter <- 1000
  
  # when there is a range of chromosomes sample a single chromosome
  dat <- chromSampler(dat.initial)
  
  
  
  # trait and data overlap
  new.dat <- TraitOverlap(dat = dat,
                          trees = trees,
                          hosts = hosts)
  
  # this dataset has diploid chromosome numbers. make them haploid
  new.dat$haploid <- new.dat$haploid / 2 
  
  # rename tree tips to match with trait data set
  tree <- trees[[i]]
  tree <- changeTipnames(tree = tree, dat = new.dat)
  
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
  chars <- new.dat$hosts
  names(chars) <- new.dat$species
  
  
  simmaps[[i]] <- make.simmap(tree = tree,
                              x = chars,
                              pi = "estimated",
                              nsim = 100,
                              message = TRUE,
                              model =  "ARD")
  
  simmap.summaries[[i]] <- describe.simmap(simmaps[[i]])
  
}

