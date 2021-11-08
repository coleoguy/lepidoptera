# Terrence Sylvester
# pradakshanas@gmail.com
# 24 June 2020

# load libraries
library(phytools)

# read data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)
# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

#seperating family genus and species names
sep.hosts <- gregexpr(pattern = "_",
                text = hosts$sp)

# get the entries of species with no family names
no.family.data.for.hosts <- c()

for(i in 1:length(sep.hosts)){
  if(length(attr(sep.hosts[[i]], which = "match.length")) < 2) no.family.data.for.hosts <- c(no.family.data.for.hosts, i)
}

# make a new matrix
dat <- matrix(data = NA,
              nrow = nrow(hosts),
              ncol = 4)

colnames(dat) <- c("family",
                   "genus",
                   "species",
                   "type")

for(i in 1:nrow(hosts)){
  if(i %in% no.family.data.for.hosts){
   entry.name <- unlist(strsplit(hosts$sp[i], split = "_"))
   
   dat[i,1] <- NA
   dat[i,2] <- entry.name[1]
   dat[i,3] <- entry.name[2]
   dat[i,4] <- hosts$hosts[i]
   
  }else{
    entry.name <- unlist(strsplit(hosts$sp[i], split = "_"))
    
    dat[i,1] <- entry.name[1]
    dat[i,2] <- entry.name[2]
    dat[i,3] <- entry.name[3]
    dat[i,4] <- hosts$hosts[i]
  }
}

# save file
write.csv(x = dat,
          file = "../data/hosts/papilionoidea-hosts.csv",
          row.names = F)

# clear envisonment
rm(list = ls())

# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

#seperating family genus and species names

for(i in 1:length(trees)){
  sep.tree <- gregexpr(pattern = "_",
                       text = trees[[i]]$tip.label)
  
  # get the entries of species with no family names
  no.family.data.for.trees <- c()
  
  for(j in 1:length(sep.tree)){
    if(length(attr(sep.tree[[j]], which = "match.length")) < 2) {
      trees[[i]]$tip.label[j] <-  paste(unlist(strsplit(x = trees[[i]]$tip.label[j], split = "_"))[1],
                                        unlist(strsplit(x = trees[[i]]$tip.label[j], split = "_"))[2],
                                        sep = " ")
    }else{
      trees[[i]]$tip.label[j] <- paste(unlist(strsplit(x = trees[[i]]$tip.label[j], split = "_"))[2], 
                                       unlist(strsplit(x = trees[[i]]$tip.label[j], split = "_"))[3], 
                                       sep = " ")
    }
  }
}

write.tree(phy = trees, 
           file = "../data/trees/processed.trees.new")
