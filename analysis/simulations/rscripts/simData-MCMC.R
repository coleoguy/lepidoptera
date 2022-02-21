# Terrence Sylvester
# pradakshanas@gmail.com
# January 25th 2022

# load libraries
library(chromePlus)
library(diversitree)
library(doSNOW)

# get helper functions
source("functions.R")
# load data
load("../results/simData/simData.RData")
# define number of clusters for parallel computing
NumberOfClusters <- 100
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)
# run MCMC
for(i in 1:length(simDat)){                      # access conditions
  for(ii in 1:length(simDat[[1]][[1]])){         # access ntips
    for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
      # make an object to store results
      results <- NULL
      results <- vector(mode = "list", length = 100)
      file.name <- NULL
      cond <- ntips <- rr <-  NULL
      cond <- names(simDat)[i]
      ntips <- names(simDat[[cond]][[1]])[ii]
      rr <- names(simDat[[cond]][[1]][[ntips]])[iii]
      file.name <- paste(cond,
                         ntips,
                         rr,
                         sep = ".")
      results <- foreach(iiii = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
        print(paste("tree", iiii))
        x <-  runMCMC(tree = trees[[ii]][[iiii]],
                      chroms = simDat[[i]][[1]][[ii]][[iii]][[iiii]],
                      binary = simDat[[i]][[2]][[ii]][[iii]][[iiii]],
                      args.lik = list(control = "ode",
                                      strict = F),
                      args.conlik = list(hyper = T,
                                         polyploidy = F,
                                         verbose = F,
                                         oneway = F,
                                         drop.poly=T,
                                         drop.demi=T,
                                         symmetric=F,
                                         nometa=F,
                                         meta="ARD"),
                      args.MCMC = list(iter.temp = 20,
                                       iter = 100,
                                       prior = make.prior.exponential(r = .5),
                                       print.every=50))
      }
      print(paste(file.name, "complete"))
      if(dir.exists("../results/simData-MCMC") == F){
        dir.create("../results/simData-MCMC")
      }
      save(results, file = paste("../results/simData-MCMC/",
                                 file.name,
                                 ".RData",
                                 sep = ""))
    }
  }
}
# stop the cluster
stopCluster(cl)
