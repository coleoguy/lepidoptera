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

anova.results <- vector(mode = "list", length = 10)

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
  
  lik <- make.mkn(tree = tree,
                  states = chrom.mat,
                  k = ncol(chrom.mat),
                  strict = FALSE,
                  control = list(method = "ode"))
  
  # model with polyploidy
  con.likwp <- constrainMkn(data = chrom.mat, 
                            lik = lik, 
                            polyploidy = F, 
                            hyper = T,
                            constrain=list(drop.demi=T))
  argnames(con.likwp)
  
  
  # model without polyploidy
  con.likwop <- constrainMkn(data = chrom.mat, 
                             lik = lik, 
                             polyploidy = F, 
                             hyper = T,
                             constrain=list(drop.demi=T, 
                                            drop.poly=T))
  
  argnames(con.likwop)
  
  fitw <- find.mle(con.likwp, 
                   x.init=c(rep(0.1, length(argnames(con.likwp)))), 
                   method = "subplex", 
                   lower = 0, 
                   upper = 145)
  
  fitwo <- find.mle(con.likwop,
                    x.init=c(rep(0.1, length(argnames(con.likwop)))),
                    method = "subplex",
                    lower = 0, 
                    upper = 145)  
  
  anova.results[[i]] <- anova(fitw,fitwo)

  
}

save.image("../results/likelihood.ratio.test.results.Rdata")

# note: only four trees were successfull. others gave an error massege.
#       based on the four trees full model is preferred