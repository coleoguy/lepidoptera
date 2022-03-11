# Terrence Sylvester
# pradakshanas@gmail.com
# january 25th 2022

# load libraries
library(chromePlus)
library(diversitree)

# get helper functions
source("functions.R")

# Build 100 trees of sizes 50,100, and 200 species --------------------
taxa <- c(50, 100, 200)
trees <- simTrees(ntaxa = taxa,
                  ntrees = 100,
                  birth = 3,
                  death = 1)
# Build 100 trees of sizes 50,100, and 200 species --------------------

# Simulate chromosome and binary trait data ---------------------
# 100 datasets for every situation
# 1) signal with binary trait no fast tips
# 2) signal with binary trait single fast tip
# 3) signal with binary trait 10% fast tips
# set rates ratio
rr <- c(1,2,5,10)
# make a place holder for simulated data
rr1 <- vector(mode = "list", length = 100)
names(rr1) <- paste("tree", 1:100, sep = "")
rr2 <- rr5 <- rr10 <- rr1
# make a place holder for rr
nTips50 <- list(rr1 = rr1,
                rr2 = rr2,
                rr5 = rr5,
                rr10 = rr10)
nTips100 <- nTips200 <- nTips50
# make a place holder for trees
chroms <- binary <- list(nTips50 = nTips50,
                         nTips100 = nTips100,
                         nTips200 = nTips200)
cond.1 <- list(chroms = chroms,
               binary = binary)
cond.2 <- cond.3 <- cond.1
# make the final dataset
simDat <- list(cond.1 = cond.1,
               cond.2 = cond.2,
               cond.3 = cond.3)
# remove unwanted data
rm(rr1,rr2,rr5,rr10,
   nTips50,nTips100,nTips200,
   chroms,binary,
   cond.1,cond.2,cond.3)

# set different rates for chromosome number evolution
chromRate <- c(0.25,0.5,0.75,1,1.25,1.5,1.75,2, 2.25,2.5,2.75,3)

ChromRanges <- vector(mode = "list", length = length(chromRate))
names(ChromRanges) <- chromRate
rngCounter <- 1

for(k in chromRate){
  
  # cycle through tree sizes
  for(i in 1:length(simDat)){                      # access conditions
    for(ii in 1:length(simDat[[1]][[1]])){         # access ntips
      for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
        for(iiii in 1:100){                        # access each tree
          check <- F
          if(i == 1){ # condition 1
            # get the tree for simulating data
            sim.tree <- NULL
            sim.tree <- trees[[ii]][[iiii]]
          }
          if(i == 2){ # condition 2
            # randomly sample a tip 
            rTip <- NULL
            rTip <- sample(Ntip(trees[[ii]][[iiii]]),1)
            # get the tree for simulating data
            sim.tree <- NULL
            sim.tree <- trees[[ii]][[iiii]]
            # increase the sample tip branch by a unit length
            sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] + 1
          }
          if(i == 3){ # condition 3
            # randomly sample a tip 
            rTip <- NULL
            rTip <- sample(Ntip(trees[[ii]][[iiii]]),
                           Ntip(trees[[ii]][[iiii]]) * 0.1,
                           replace = F)
            # get the tree for simulating data
            sim.tree <- NULL
            sim.tree <- trees[[ii]][[iiii]]
            # increase the sampled tip branches by a unit length
            sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] + 1
          }
          # simulate chromosomes and binary states
          while(check==F){
            z <- NULL
            z <- simChrom(sim.tree, 
                          pars=c(k,         # gains at state 0
                                 rr[iii]*k, # gains at state 1
                                 k,         # loss at state 0
                                 rr[iii]*k, # loss at state 1
                                 0,           # demiploidy at state 0
                                 0,           # demiploidy at state 1
                                 0,           # polyploidy at state 0
                                 0,           # polyploidy at state 1
                                 .5,          # transition from state 0 to state 1
                                 .5,          # transition from state 1 to state 0
                                 10,          # root chromosome number
                                 0),          # root state
                          limits = c(1, 100), model = "ChromPlus")
            if(sum(z$binary.state) >= (taxa[ii]*0.10) &
               sum(z$binary.state) <= (taxa[ii]*0.90)){
              check <- T
            }
          }
          # store chroms
          simDat[[i]][[1]][[ii]][[iii]][[iiii]] <- z[[2]]
          # store binary trait
          simDat[[i]][[2]][[ii]][[iii]][[iiii]] <- z[[1]]
        }
      }
    }
  }
  
  # make a table to store chrom range data
  range.stat <-  as.data.frame(matrix(data = NA,
                                      nrow = 36,
                                      ncol = 4))
  colnames(range.stat) <- c("ntips", "condition", "rr", "chromRange")
  # add a counter
  counter <- 1
  # cycle through simDat to find median range
  for(i in 1:length(simDat)){                      # access conditions
    for(ii in 1:length(simDat[[1]][[1]])){         # access ntips
      for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
        rng <- NULL
        rng <- vector(mode = "numeric", length = 100)
        for(iiii in 1:100){
          rng.1 <- range(simDat[[i]][[1]][[ii]][[iii]][[iiii]])[1]
          rng.2 <- range(simDat[[i]][[1]][[ii]][[iii]][[iiii]])[2]
          rng[iiii] <- rng.2 - rng.1
        }
        range.stat$condition[counter] <- as.numeric(gsub(pattern = "cond.", replacement = "", x = names(simDat)[i]))
        range.stat$ntips[counter] <- as.numeric(gsub(pattern = "nTips", replacement = "", x = names(simDat[[1]][[1]])[ii]))
        range.stat$rr[counter] <- as.numeric(gsub(pattern = "rr", replacement = "", x = names(simDat[[1]][[1]][[1]])[iii]))
        range.stat$chromRange[counter] <- median(rng)
        counter <- counter + 1
      }
    }
  }
  range.stat[order(range.stat$ntips),]
  ChromRanges[[rngCounter]] <- range.stat
  rngCounter <- rngCounter + 1
}

par(mfrow = c(3,4))
for(ii in 1:length(simDat)){                      # access conditions
  for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
    
    cond <- as.numeric(gsub(pattern = "cond.", replacement = "", x = names(simDat)[ii]))
    rr <- as.numeric(gsub(pattern = "rr", replacement = "", x = names(simDat[[1]][[1]][[1]])[iii]))
    
    # par(xaxt = "n")
    plot(x = NULL,
         y = NULL,
         xlim = c(0,3),
         ylim = c(0,20),
         xlab = "Base rate",
         ylab = "Median range",)
    # par(xaxt = "s")
    # axis(side = 1,at = names(ChromRanges),
    #      labels = names(ChromRanges),
    #      tick = T,srt = -90)
    # # text(names(ChromRanges),
    #      par("usr")[3],
    #      labels = names(ChromRanges),
    #      srt = -90, 
    #      pos = 1,
    #      xpd = T,
    #      adj = 0.5,
    #      offset = 1,)
    
    
    mtext(text = paste("Condition =", cond),side = 3,adj = 0,
          cex = 0.7)
    mtext(text = paste("Rates ratio =", rr),side = 3,adj = 1,
          cex = 0.7)
    
    abline(h = 8, col = "gray", lty = 2)
    
    for(i in 1:length(ChromRanges)){
      points(x = as.numeric(names(ChromRanges)[i]),
             y = ChromRanges[[i]]$chromRange[ChromRanges[[i]]$ntips == 50 & ChromRanges[[i]]$condition == cond & ChromRanges[[i]]$rr == rr],
             pch = 16,
             col = rgb(1,0,0,0.5))  
    }
    
    for(i in 1:length(ChromRanges)){
      points(x = as.numeric(names(ChromRanges)[i]),
             y = ChromRanges[[i]]$chromRange[ChromRanges[[i]]$ntips == 100 & ChromRanges[[i]]$condition == cond & ChromRanges[[i]]$rr == rr],
             pch = 16,
             col = rgb(0,1,0,0.5))  
    }
    
    for(i in 1:length(ChromRanges)){
      points(x = as.numeric(names(ChromRanges)[i]),
             y = ChromRanges[[i]]$chromRange[ChromRanges[[i]]$ntips == 200 & ChromRanges[[i]]$condition == cond & ChromRanges[[i]]$rr == rr],
             pch = 16,
             col = rgb(0,0,1,0.5))  
    }
  }
}

legend("bottomright",inset = 0.02,
       legend = c(50,100,200),
       pch = 16,
       col = c("red", "green","blue"),)

