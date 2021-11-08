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

simmap.otto <- vector(mode = "list", length = 10)
names(simmap.otto) <- paste("tree", 1:10, sep = ".")


# make q matrix using the rates from Hardy and Otto 2014 paper

## speciation and extinction rates are different in specialist vs generalist species
qmat <- matrix(data = c(-0.13726719, 0.13726719, 0.08788453, -0.08788453),
               nrow = 2,
               ncol = 2,
               byrow = T,dimnames = list(c(0,1),c(0,1)))

for(i in 1:10){
simmap.otto[[i]] <- make.simmap(tree = trees[[i]],
            x = feedingType,
            Q = qmat,
            pi = c(1,0),
            nsim = 100)
}

save.image("../results/stoch.mapping.otto.RData")
