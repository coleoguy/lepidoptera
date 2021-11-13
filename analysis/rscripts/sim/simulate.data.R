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
rm(tree, tree.set, i, j, taxa)
# Build 100 trees of sizes 50 to 250 species --------------------


# Simulate chromosome and binary trait data ---------------------
# 100 datasets for every situation
# 1) no signal with binary trait no violation
# 2) signal with binary trait no violation
# 3) no signal with binary trait single tip violation
# 4) signal with binary trait single tip violation
# 5) one size of tree but increase number of tip violators
#    up to 10% of tips

# for sims with signal rates ratio will be 5 to replicate 2019 paper
# for sims with single tip violation branch will be extended by 1
# unit

# this will contain the five datasets described above
chrom.traits <- bin.traits <- list()

# cycle through tree sizes
for(i in 1:5){
  # condition 1
  for(j in 1:100){
    check <- F
    while(check==F){
      total.counter <- total.counter + 1
      z <- simChrom(trees[[i]][[j]], 
                    pars=c(.75, 1*.75, .75, 1*.75,
                           0,   0,  0,  0, .5,  .5, 10,  0),
                    limits = c(1, 100), model = "ChromPlus")
      if(sum(z$binary.state) > 5 &
         sum(z$binary.state) < (length(z$binary.state)-5)){
        check <- T
      }
    }
    bin.traits

  }
  
  # condition 2
  for(j in 1:100){
    
  }
  
  
  # condition 3
  for(j in 1:100){
    
  }
  
  
  # condition 4
  for(j in 1:100){
    
  }
  
  
  # condition 5
  for(j in 1:100){
    
  }
  
  
  
  
  
  
  
  
}




# Simulate chromosome and binary trait data ---------------------




