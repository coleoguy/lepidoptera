# Terrence Sylvester
# 13 October 202
# pradakshanas@g,ail.com

# load libraries
library(phytools)

# read in data
trees <- read.tree("../data/trees/processed.trees.new")
hosts <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)

# combine genus and species names in hosts
hosts$spNames <- paste(hosts$genus, hosts$species, sep = "_")

# feeding type
feedingType <- hosts$type
names(feedingType) <- hosts$spNames

simmap <- vector(mode = "list", length = 10)
names(simmap) <- paste("tree", 1:10, sep = ".")

for(i in 1:10){
simmap[[i]] <- make.simmap(tree = trees[[i]],
            x = feedingType,
            model = "ARD",
            pi = c(1,0),
            nsim = 100)
}

save.image("../results/stoch.mapping.RData")
