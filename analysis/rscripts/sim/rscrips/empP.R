# load libraries
library(phytools)
library(diversitree)
library(chromePlus)

# load helper functions
source("functions.R")
# load data
load("simData.RData")
load("cond.1.nTips50.rr1.RData")
# get post burnin of emperical data
pbrn <- getPostBurnin(results, burn = 0.5)
# run MCMC
res <- getEmpiricalP(tree = trees$nTips50[[1]],
           chroms = simDat$cond.1$chroms$nTips50$rr1$tree1,
           binary = simDat$cond.1$binary$nTips50$rr1$tree1,
           data = pbrn[1:50,],
           nsim = 1,
           plot.lik = T,
           plot.p = T)