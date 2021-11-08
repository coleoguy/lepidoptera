# Terrence Sylvester
# pradakshanas@gmail.com
# 2nd Oct, 2019

# load libraries
library(diversitree)

# read in data
chroms <- read.csv("../data/chroms/lep.chroms.csv", as.is = T)

# read in trees
trees <- read.tree('../Hardy-Otto-Data 2/Papilionoidea-BiSSE/papilionoidea-10.trees')

# read in hosts data
d <- read.csv('../Hardy-Otto-Data 2/Papilionoidea-BiSSE/papilionoidea-bisse-data.txt')

# compare the tree datasets
lep.names <- chroms$species
tree.names <- trees[[1]]$tip.label
d.names <-  as.character(d$sp)

# get the overlap between the tree datasets
sum(unique(lep.names) %in%  unique(tree.names))

sum(unique(tree.names) %in% unique(lep.names))

sum(unique(tree.names) %in% unique(d.names))

# number of specialist species
sum(d$hosts[d.names %in% unique(lep.names)])