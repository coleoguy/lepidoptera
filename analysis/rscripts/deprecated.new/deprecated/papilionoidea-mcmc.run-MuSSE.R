# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# mcmc analysus

# load libraries
library(diversitree)
library(chromePlus)
library(phytools)

# load helper functions
source("helper.functions.R")
# read in chromosome data
dat.initial <- read.csv("../data/chroms/lep.chroms.csv", as.is = T)
# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")
# read in hosts data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)
results <- vector(mode = "list", length = 10)
tree.depth <- vector(mode = "numeric", length = 10)
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
  tree.depth[i] <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth[i]
  # track progress
  print(paste("working on tree ", i))
  rng <- c(range(new.dat$haploid, na.rm = T)[1] - 2,
           range(new.dat$haploid, na.rm = T)[2] + 2)
  # here as the hyper trait this should be the probability in being in state 1
  # of trait 2. Which means the probability of being a generalist. Because
  # our state 1 is generalists and state 2 is specialist. We have to change
  # our hosts dataset as this shows the probability of being generalist as
  # state 1. We change that here
  binary.trait <- new.dat$hosts == 1
  new.dat$hosts[binary.trait == T] <- 0
  new.dat$hosts[binary.trait == F] <- 1
  chrom.mat <- datatoMatrix(x = new.dat[,c(1,2,3)], 
                            range = rng,
                            hyper = T)
  lik <- make.musse(tree = tree,
                    states = chrom.mat,
                    k = ncol(chrom.mat),
                    strict = FALSE,
                    control = list(method = "ode"))
  con.lik <- constrainMuSSE(data = chrom.mat,
                            lik = lik,
                            polyploidy = F,
                            s.lambda = F,
                            s.mu = F,
                            hyper = T,
                            constrain = list(drop.demi = T,
                                             drop.poly = F))
  results[[i]] <- mcmc(lik = con.lik,
                       x.init = rep(0.1, length(argnames(con.lik))),
                       nsteps = iter,
                       w = 1,
                       prior= make.prior.exponential(r=10))
  # here asc 1 and desc 1 are fusoin and fissoin rates of generalists
  # and asc 2 and esc 2 are fusion and fussion rates of specialists
  # convert rates to millions of years
  results[[i]][,c(1:length(argnames(con.lik)))+1] <- results[[i]][,c(1:length(argnames(con.lik)))+1] / tree.depth[i]
}
