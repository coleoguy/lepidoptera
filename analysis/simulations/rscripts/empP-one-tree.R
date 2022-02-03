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
# make a place holder for emperical data
empP <- vector(mode = "list", length = 100)
names(empP) <- paste("tree", 1:100, sep = ".")

# get list of completed trees
compTrees <- as.numeric(gsub(pattern = "tree.", replacement = "", dir("../results/empP-MCMC/")))
check <- T
while (check) {
  # sample a single tree
  iiii <- sample(1:100,1)
  if(iiii %in% compTrees == F)
  check <- F
}

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
      # get post burnin of emperical data
      pbrn <- getPostBurnin(results[[iiii]], burn = 0.5)
      # run MCMC
      empP <- getEmpiricalPMC(tree = trees[[ii]][[iiii]],
                              chroms = simDat[[i]][[1]][[ii]][[iii]][[iiii]],
                              binary = simDat[[i]][[2]][[ii]][[iii]][[iiii]],
                              data = pbrn,
                              nsim = 100,
                              plot.lik = T,
                              plot.p = T,
                              nclust = 100,
                              burn = 0.5)
      if(dir.exists(paste("../results/empP-MCMC/tree",iiii, sep=".")) == F){
        dir.create(paste("../results/empP-MCMC/tree",iiii, sep="."))  
      }
      save(empP, file = paste("../results/empP-MCMC/tree.",
                              iiii,
                              "/empP",
                              file.name,
                              ".tree.",
                              iiii,
                              ".RData", sep = ""))
      # remove emp from memory
      rm(empP)
      # free unused memory
      gc()
      print(paste(file.name, "complete"))
    }
    # remove results
    rm(results)
    # free unused memory
    gc()
  }
}

