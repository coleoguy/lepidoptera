# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# mcmc analysis

# load libraries
library(diversitree)
library(chromePlus)
library(phytools)
library(doMC) # this allows multicore runs on a mac

# load helper functions
source("helper.functions.R")

# read in chromosome data
dat.initial <- read.csv("../data/chroms/lep.chroms.csv", as.is = T)
# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")
# read in hosts data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)
tree.depth <- vector(mode = "numeric", length = 100)



registerDoMC(14)
####result <- list()

# we will loop through all 100 trees
# fitting model without polyploidy
iter <- 100
iter.temp <- 20
tree.num <- rep(1:10, each=10)
prior <- make.prior.exponential(r=2)
upper <- 50

# get w
dat <- chromSampler(dat.initial)
# trait and data overlap
new.dat <- TraitOverlap(dat = dat,
                        trees = trees,
                        hosts = hosts)
new.dat$haploid <- new.dat$haploid / 2
tree <- trees[[tree.num[i]]]
# this was a confusing name to give your function this
# where you are pruning your trees too. Most people
# would expect that to be in trait overlap
tree <- changeTipnames(tree = tree, dat = new.dat)
tree <- force.ultrametric(tree = tree)
rng <- c(range(new.dat$haploid, na.rm = T)[1] - 1,
         range(new.dat$haploid, na.rm = T)[2] + 1)
# our states are coded such that specialists will be state 1 and
# generalists will be state 2
chrom.mat <- datatoMatrix(x = new.dat[,1:3],
                          range = rng,
                          hyper = T)
# make the likelihood function
lik <- make.mkn(tree = tree,
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = FALSE,
                control = list(method = "ode"))
# constrain the likelihood function
con.lik <- constrainMkn(data = chrom.mat,
                        lik = lik,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly = F))
temp <- mcmc(lik = con.lik,
             x.init = runif(min=0, max=1,
                            n=length(argnames(con.lik))),
             nsteps = iter.temp,
             w = 1,
             prior = prior,
             upper = upper)
w <- diff(sapply(temp[11:20, 2:9], quantile, c(.05, .95)))




x <- foreach(i=1:100) %dopar%{
  # when there is a range of chromosomes sample a single chromosome
  dat <- chromSampler(dat.initial)
  # trait and data overlap
  new.dat <- TraitOverlap(dat = dat,
                          trees = trees,
                          hosts = hosts)
  # this dataset has diploid chromosome numbers. make them haploid
  new.dat$haploid <- new.dat$haploid / 2
  # rename tree tips to match with trait data set
  tree <- trees[[tree.num[i]]]
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
  # print(paste("working on tree ", i))
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
  chrom.mat <- datatoMatrix(x = new.dat[,1:3],
                            range = rng,
                            hyper = T)
  # make the likelihood function
  lik <- make.mkn(tree = tree,
                  states = chrom.mat,
                  k = ncol(chrom.mat),
                  strict = FALSE,
                  control = list(method = "ode"))
  # constrain the likelihood function
  con.lik <- constrainMkn(data = chrom.mat,
                          lik = lik,
                          polyploidy = F,
                          hyper = T,
                          constrain = list(drop.demi = T,
                                           drop.poly = F))
  # run the mcmc chain
  results <- mcmc(lik = con.lik,
                  x.init = runif(min=0, max=1,
                                 n=length(argnames(con.lik))),
                  nsteps = iter,
                  w = w, #
                  prior = prior,
                  upper = upper)
  # here asc 1 and desc 1 are fusion and fission rates of generalists
  # and asc 2 and desc 2 are fusion and fission rates of specialists
  # convert rates to millions of years
  results
}


for(i in 1:100){
  dat <- chromSampler(dat.initial)
  # trait and data overlap
  new.dat <- TraitOverlap(dat = dat,
                          trees = trees,
                          hosts = hosts)
  tree <- trees[[tree.num[i]]]
  tree <- changeTipnames(tree = tree, dat = new.dat)
  tree <- force.ultrametric(tree = tree)
  # scale tree to unit length
  tree.depth[i] <- max(branching.times(tree))
}

rm(list=ls()[-c(14,18)]
# TODO
# save tree depths and raw mcmc as a result file.
# then make second script process MCMC that looks at
# burnin transforms back to MY rates and sample post burnin
