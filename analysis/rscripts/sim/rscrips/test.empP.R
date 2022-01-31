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
# get the transistion matrix
tmat <- matrix(data = 0, nrow = 2, ncol = 2)
colnames(tmat) <- rownames(tmat) <- c(0,1)
# fill tmat
tmat[1,2] <- mean(pbrn$tran12[1:50])
tmat[2,1] <- mean(pbrn$tran21[1:50])
diag(tmat) <- -rowSums(tmat)
# run MCMC
res <- simBinMCMC(tree = trees$nTips50[[1]],
           chroms = simDat$cond.1$chroms$nTips50$rr1$tree1,
           binary = simDat$cond.1$binary$nTips50$rr1$tree1,tmat = tmat,
           nsim = 10,
           hyper = T,
           polyploidy = F,
           verbose = F,
           oneway = F,
           drop.poly = T,
           drop.demi = T,
           iter.temp = 20,
           iter = 100)
#Plot MCMC
plotlikMCMC(res,0.5)
# get post burnin of MCMC
simpb <-  getPostBurnin(res, 0.5)
# calculate emperical P
empiricalPcalc(empPostburnin = pbrn[1:50,],
               simPostburnin = simpb,
               polyploidy = F,
               nsim = 10,
               plot = T)
