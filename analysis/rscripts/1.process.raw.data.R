# Terrence Sylvester
# pradakshanas@gmail.com
# 24 June 2020

# load libraries
library(phytools)

# read data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)
# read in trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

# separate family genus and species names from hosts database
# in the hosts database the species are labelled as follows: family_genus_species
# Examining the dataset before processing it shows that there are no subspecies
# information in this dataset
get.higher.taxa <- gregexpr(pattern = "_",
                            text = hosts$sp)

# that can me confirmed by forllowing script
for(i in 1:length(get.higher.taxa)){
  if(length(attr(get.higher.taxa[[i]], which = "match.length")) > 2){
    print(i)
  } 
}

# in some instances for some species the family information is left out
# get the entries of species with no family names
miss.family <- c()
for(i in 1:length(get.higher.taxa)){
  if(length(attr(get.higher.taxa[[i]], which = "match.length")) < 2) {
    miss.family <- c(miss.family, i)
  }
}

# test if this is correct
hosts$sp[miss.family]

# make a new matrix to create a new dataset
dat <- as.data.frame(matrix(data = NA,
                            nrow = nrow(hosts),
                            ncol = 5))
# get the column names
colnames(dat) <- c("family",
                   "genus",
                   "species",
                   "binomial",
                   "type")
# fill in the new dataset
for(i in 1:nrow(hosts)){
  if(i %in% miss.family){
    entry.name <- unlist(strsplit(hosts$sp[i], split = "_"))
    
    dat$family[i] <- NA
    dat$genus[i] <- entry.name[1]
    dat$species[i] <- entry.name[2]
    dat$type[i] <- hosts$hosts[i]
    
  }else{
    entry.name <- unlist(strsplit(hosts$sp[i], split = "_"))
    
    dat$family[i] <- entry.name[1]
    dat$genus[i] <- entry.name[2]
    dat$species[i] <- entry.name[3]
    dat$type[i] <- hosts$hosts[i]
    }
}
# fill in the binomial
# make sure that the seperator between genus and species names matches with that
# of the phylogeny
trees[[1]]$tip.label[1]
for(i in 1:nrow(dat)){
  dat$binomial[i] <- paste(dat$genus[i], dat$species[i], sep = "_")
}

#separating family genus and species names
for(i in 1:length(trees)){
  sep.tree <- gregexpr(pattern = "_",
                       text = trees[[i]]$tip.label)
  # remove the family level information from the tree tip lables
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

# trees are not ultrametric
# convert the trees to ultrametric form
phy <- vector(mode = "list", length = length(trees))
for(i in 1:length(trees)){
  phy[[i]] <- force.ultrametric(trees[[i]])
}
# change the class of the new set of trees to match with the original set of trees
class(phy) <- class(trees)

# save hosts and trees to a file
write.csv(x = dat,
          file = "../data/hosts/papilionoidea-hosts.csv",
          row.names = F)
write.tree(phy = phy,
           file = "../data/trees/processed.trees.new")
