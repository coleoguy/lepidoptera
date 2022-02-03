# Terrence Sylvester
# pradakshanas@gmail.com
# january 25th 2022

# load libraries
library(chromePlus)
library(diversitree)

# set initial conditions
# set rates ratio
rr <- c(1,2,5,10)
# set numer of taxa
taxa <- c(50, 100, 200)

# Build 100 trees of sizes 50,100, and 200 species --------------------
trees <- list()

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
names(trees) <- paste("nTips",taxa, sep = "") 
# remove unwanted objects
rm(tree, tree.set, i, j)
# Build 100 trees of sizes 50,100, and 200 species --------------------

# Simulate chromosome and binary trait data ---------------------
# 100 datasets for every situation

# 1) signal with binary trait no fast tips
# 2) signal with binary trait single fast tip
# 3) signal with binary trait 10% fast tips

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

# cycle through tree sizes
for(i in 1:length(simDat)){                      # access conditions
  for(ii in 1:length(simDat[[1]][[1]])){         # access ntips
    for(iii in 1:length(simDat[[1]][[1]][[1]])){ # access rr
      for(iiii in 1:100){                        # access each tree
        check <- F
        if(i == 1){ # simulate data for condition 1
          while(check==F){
            z <- NULL
            z <- simChrom(trees[[ii]][[iiii]], 
                          pars=c(.75,         # gains at state 0
                                 rr[iii]*.75, # gains at state 1
                                 .75,         # loss at state 0
                                 rr[iii]*.75, # loss at state 1
                                 0,       # demiploidy at state 0
                                 0,       # demiploidy at state 1
                                 0,       # polyploidy at state 0
                                 0,       # polyploidy at state 1
                                 .5,      # transition from state 0 to state 1
                                 .5,      # transition from state 1 to state 0
                                 10,      # root chromosome number
                                 0),      # root state
                          limits = c(1, 100), model = "ChromPlus")
            if(sum(z$binary.state) >= (taxa[ii]*0.10) &
               sum(z$binary.state) <= (taxa[ii]*0.90)){
              check <- T
            }
          }
        }
        if(i == 2){ # simulate data for condition 2
          # randomly sample a tip 
          rTip <- NULL
          rTip <- sample(Ntip(trees[[ii]][[iiii]]),1)
          # get the tree for simulating data
          sim.tree <- NULL
          sim.tree <- trees[[ii]][[iiii]]
          # increase the sample tip branch by a unit length
          sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] == rTip)] + 1
          while(check==F){
            z <- NULL
            z <- simChrom(sim.tree, 
                          pars=c(.75,     # gains at state 0
                                 rr[iii]*.75,   # gains at state 1
                                 .75,     # loss at state 0
                                 rr[iii]*.75,   # loss at state 1
                                 0,       # demiploidy at state 0
                                 0,       # demiploidy at state 1
                                 0,       # polyploidy at state 0
                                 0,       # polyploidy at state 1
                                 .5,      # transition from state 0 to state 1
                                 .5,      # transition from state 1 to state 0
                                 10,      # root chromosome number
                                 0),      # root state
                          limits = c(1, 100), model = "ChromPlus")
            if(sum(z$binary.state) >= (taxa[ii]*0.10) &
               sum(z$binary.state) <= (taxa[ii]*0.90)){
              check <- T
            }
          }
        }
        if(i == 3){ # simulate data for condition 3
          # randomly sample a tip 
          rTip <- NULL
          rTip <- sample(Ntip(trees[[ii]][[iiii]]),10, replace = F)
          # get the tree for simulating data
          sim.tree <- NULL
          sim.tree <- trees[[ii]][[iiii]]
          # increase the sampled tip branches by a unit length
          sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] <- sim.tree$edge.length[which(sim.tree$edge[,2] %in% rTip)] + 1
          while(check==F){
            z <- NULL
            z <- simChrom(sim.tree, 
                          pars=c(.75,     # gains at state 0
                                 rr[iii]*.75,   # gains at state 1
                                 .75,     # loss at state 0
                                 rr[iii]*.75,   # loss at state 1
                                 0,       # demiploidy at state 0
                                 0,       # demiploidy at state 1
                                 0,       # polyploidy at state 0
                                 0,       # polyploidy at state 1
                                 .5,      # transition from state 0 to state 1
                                 .5,      # transition from state 1 to state 0
                                 10,      # root chromosome number
                                 0),      # root state
                          limits = c(1, 100), model = "ChromPlus")
            if(sum(z$binary.state) >= (taxa[ii]*0.10) &
               sum(z$binary.state) <= (taxa[ii]*0.90)){
              check <- T
            }
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

# remove unwanted objects from the environment
rm(i,ii,iii,iiii,z, check, rTip, sim.tree)

# save data
save(simDat,trees, file = "simData.RData")