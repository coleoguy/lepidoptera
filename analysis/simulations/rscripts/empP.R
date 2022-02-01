# load libraries
library(phytools)
library(diversitree)
library(chromePlus)

# load helper functions
source("functions.R")
# load data
load("simData.RData")
# make a place holder for emperical data
empP <- vector(mode = "list", length = 100)
names(empP) <- paste("tree", 1:100, sep = ".")
# get results of each condition
for(i in 1:length(simDat)){ # access conditions
  for(ii in 1:length(simDat[[1]][[1]])){ # access ntips
    for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
      cond <- names(simDat)[i]
      ntips <- names(simDat[[1]][[1]])[ii]
      rr <- names(simDat[[1]][[1]][[1]])
      # get the save file name
      file.name <- paste(cond,ntips,rr,"RData",sep = ".")
      #load results
      load(file.name)
      for(iiii in 1:length(results)){
        # get post burnin of emperical data
        pbrn <- getPostBurnin(results[[iiii]], burn = 0.5)
        # run MCMC
        empP[[iiii]] <- getEmpiricalP(tree = trees[[ii]][[iiii]],
                             chroms = simDat[[i]][[1]][[ii]][[iii]][[iiii]],
                             binary = simDat[[i]][[2]][[ii]][[iii]][[iiii]],
                             data = pbrn,
                             nsim = 100,
                             plot.lik = F,
                             plot.p = F)
      }
      print(paste(file.name, "complete"))
      save(empP, file = paste("empP",file.name,"tree",iiii,"RData", sep = "."))
    }
  }
}
    
    
    
    
    
    