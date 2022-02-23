# Terrence Sylvester
# pradakshanas@gmail.com
# 22nd February 2022

# load libraries
library(phytools)
library(diversitree)
library(chromePlus)
library(maps)
library(doSNOW)

# load helper functions
source("functions.R")
# load data
load("../results/simData/simData.RData")
# make a place holder for empirical data
empP <- vector(mode = "list", length = 100)
names(empP) <- paste("tree", 1:100, sep = ".")
# set pars
# get results of each condition
for(i in 1:length(simDat)){ # access conditions
  for(ii in 1:length(simDat[[1]][[1]])){ # access ntips
    for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
      cond <- names(simDat)[i]
      ntips <- names(simDat[[1]][[1]])[ii]
      rr <- names(simDat[[1]][[1]][[1]])[iii]
      # get the save file name
      file.name <- paste(cond,ntips,rr,"RData",sep = ".")
      #load results
      load(paste("../results/simData-MCMC/",file.name, sep = ""))
      for(iiii in 1:100){
        # get post burn-in of empirical data
        pbrn <- getPostBurnin(results[[iiii]], burn = 0.5)
        # run MCMC
        empP <- getEmpiricalPMC(tree = trees[[ii]][[iiii]],
                                chroms = simDat[[i]][[1]][[ii]][[iii]][[iiii]],
                                binary = simDat[[i]][[2]][[ii]][[iii]][[iiii]],
                                data = pbrn,
                                nsim = 100,
                                plot.lik = F,
                                plot.p = F,
                                nclust = 100,
                                burn = 0.5,
                                args.MCMC = list(iter.temp = 20,
                                                 iter = 100,
                                                 prior = make.prior.exponential(r = .5),
                                                 print.every=100))
        if(dir.exists(paste("../results/empP-MCMC-all-trees/tree",iiii, sep=".")) == F){
          dir.create(paste("../results/empP-MCMC-all-trees/tree",iiii, sep="."),recursive = T)  
        }
        save(empP, file = paste("../results/empP-MCMC-all-trees/tree.",
                                iiii,
                                "/empP.",
                                file.name,
                                sep = ""))
        # remove emp from memory
        rm(empP)
        # free unused memory
        gc()
      }
      print(paste(file.name, "complete"))
    }
    # remove results
    rm(results)
    # free unused memory
    gc()
  }
}
