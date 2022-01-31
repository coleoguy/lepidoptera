# This script simulates trees and chromosome number data
# It includes datasets with and without significant outlier
# taxa

library(diversitree)
library(chromePlus)


# Build 100 trees of sizes 50 to 250 species --------------------
trees <- list()
taxa <- c(50, 100, 150, 200, 250)
for(i in 1:length(taxa)){
  tree.set <- list()
  for(j in 1:100){
    print(j)
    tree <- NULL
    while(is.null(tree)){
      tree <- tree.bd(pars = c(3, 1), max.taxa = taxa[i])
    }
    tree$edge.length <- tree$edge.length/max(branching.times(tree))
    tree.set[[j]] <- tree
  }
  trees[[i]] <- tree.set 
}
# name each set in trees
names(trees) <- taxa
# remove unwanted objects
rm(tree, tree.set, i, j)
# Build 100 trees of sizes 50 to 250 species --------------------


# Simulate chromosome and binary trait data ---------------------
# 100 datasets for every situation

# 1) no signal with binary trait no violation
# 2) signal with binary trait no violation
# 3) no signal with binary trait single tip violation
# 4) signal with binary trait single tip violation
# 5) one size of tree but increase number of tip violators
#    up to 10% of tips
# 6) Signal with binary trait and one size of tree but increase number of tip violators

# for sims with signal rates ratio will be 5 to replicate 2019 paper
# for sims with single tip violation branch will be extended by 1
# unit

# this will contain the five datasets described above
chrom.traits.1 <- bin.traits.1 <- vector(mode = "list", length = length(trees))
chrom.traits.2 <- bin.traits.2 <- vector(mode = "list", length = length(trees))
chrom.traits.3 <- bin.traits.3 <- vector(mode = "list", length = length(trees))
chrom.traits.4 <- bin.traits.4 <- vector(mode = "list", length = length(trees))
chrom.traits.5 <- bin.traits.5 <- vector(mode = "list", length = 5)
chrom.traits.6 <- bin.traits.6 <- vector(mode = "list", length = 5)

# assign names for each dataset
names(chrom.traits.1) <- names(bin.traits.1) <- names(trees)
names(chrom.traits.2) <- names(bin.traits.2) <- names(trees)
names(chrom.traits.3) <- names(bin.traits.3) <- names(trees)
names(chrom.traits.4) <- names(bin.traits.4) <- names(trees)
names(chrom.traits.5) <- names(bin.traits.5) <- names(trees)
names(chrom.traits.5) <- names(bin.traits.5) <-  c(2,4,6,8,10)
names(chrom.traits.6) <- names(bin.traits.6) <-  c(2,4,6,8,10)

# set the rates ratio
rr <- 10

# cycle through tree sizes
for(i in 1:5){
  # condition 1
  # 1) no signal with binary trait no violation
  for(j in 1:100){
    check <- F
    while(check==F){
      z <- NULL
      z <- simChrom(trees[[i]][[j]], 
                    pars=c(.75,     # gains at state 0
                           1*.75,   # gains at state 1
                           .75,     # loss at state 0
                           1*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.1[[i]][[j]] <- z[[1]]
    chrom.traits.1[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 1 complete"))
  # condition 2
  # 2) signal with binary trait no violation
  for(j in 1:100){
    check <- F
    while(check==F){
      z <- NULL
      z <- simChrom(trees[[i]][[j]], 
                    pars=c(.75,     # gains at state 0
                           rr*.75,   # gains at state 1
                           .75,     # loss at state 0
                           rr*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.2[[i]][[j]] <- z[[1]]
    chrom.traits.2[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 2 complete"))
  # condition 3
  # 3) no signal with binary trait single tip violation
  
  for(j in 1:100){
    check <- F
    while(check==F){
      z <- NULL
      # randompy sample a tip 
      rTip <- NULL
      rTip <- sample(Ntip(trees[[i]][[j]]),1)
      # extend the edge of the tip
      sim.tree <- NULL
      sim.tree <- trees[[i]][[j]]
      sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] + 1
      z <- simChrom(trees[[i]][[j]], 
                    pars=c(.75,     # gains at state 0
                           1*.75,   # gains at state 1
                           .75,     # loss at state 0
                           1*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.3[[i]][[j]] <- z[[1]]
    chrom.traits.3[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 3 complete"))
  # condition 4
  # 4) signal with binary trait single tip violation
  for(j in 1:100){
    check <- F
    while(check==F){
      z <- NULL
      # randompy sample a tip 
      rTip <- NULL
      rTip <- sample(Ntip(trees[[i]][[j]]),1)
      # extend the edge of the tip
      sim.tree <- NULL
      sim.tree <- trees[[i]][[j]]
      sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] + 1
      z <- simChrom(sim.tree, 
                    pars=c(.75,     # gains at state 0
                           rr*.75,   # gains at state 1
                           .75,     # loss at state 0
                           rr*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.4[[i]][[j]] <- z[[1]]
    chrom.traits.4[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 4 complete"))
}


# condition 5
# 5) one size of tree but increase number of tip violators
#    up to 10% of tips
for(i in 1:5){
  for(j in 1:100){
    z <- NULL
    # randompy sample a tip 
    rTip <- NULL
    rTip <- sample(Ntip(trees[[5]][[j]]),Ntip(trees[[5]][[j]])*i*0.01*2,replace = F)
    # extend the edge of the tip
    sim.tree <- NULL
    sim.tree <- trees[[5]][[j]]
    sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] + 1
    check <- F
    while(check==F){
      z <- simChrom(trees[[5]][[j]], 
                    pars=c(.75,     # gains at state 0
                           1*.75,   # gains at state 1
                           .75,     # loss at state 0
                           1*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.5[[i]][[j]] <- z[[1]]
    chrom.traits.5[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 5 complete"))
}



# condition 6
# 6) Signal with Bbinary trait and one size of tree but increase number of tip violators
#    up to 10% of tips
for(i in 1:5){
  for(j in 1:100){
    z <- NULL
    # randompy sample a tip 
    rTip <- NULL
    rTip <- sample(Ntip(trees[[5]][[j]]),Ntip(trees[[5]][[j]])*i*0.01*2,replace = F)
    # extend the edge of the tip
    sim.tree <- NULL
    sim.tree <- trees[[5]][[j]]
    sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] + 1
    check <- F
    while(check==F){
      z <- simChrom(trees[[5]][[j]], 
                    pars=c(.75,     # gains at state 0
                           rr*.75,   # gains at state 1
                           .75,     # loss at state 0
                           rr*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) >= (taxa[i]*0.10) &
         sum(z$binary.state) <= (taxa[i]*0.90)){
        check <- T
      }
    }
    bin.traits.6[[i]][[j]] <- z[[1]]
    chrom.traits.6[[i]][[j]] <- z[[2]]
    
  }
  print(paste("simulations for condidion 5 complete"))
}
# Simulate chromosome and binary trait data ---------------------

# save data
save.image("simulate.data.rr-10.RData")



