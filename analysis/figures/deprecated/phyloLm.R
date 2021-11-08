# Terrence Sylvester
# pradakshanas@gmail.com
# 24 November 2020

# library
library(phytools)
library(phylolm)
# get helper functions
source("../rscripts/helper.functions.R")
# read in data
bsz <- read.csv("../data/body.size/body.size.csv", as.is = T)
bsz$speciesNames <- paste(bsz$Genus, bsz$Species, sep = "_")
# read trees
trees <- read.tree("../data/trees/processed.trees.new")
# read in chromosome number data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
# fill in species for data
dat$SpecisName <- paste(dat$other.names.genera,
                        dat$other.names.species,
                        sep = "_")
dat$SpecisName[dat$SpecisName == "_"] <- paste(dat$Genus[dat$SpecisName == "_"],
                                               dat$species[dat$SpecisName == "_"],
                                               sep = "_")
# sample a single chromosome number when there is a range
dat <- chromSampler(dat = dat)
# remove species that lack male chromosome number data
dat <-  dat[!(is.na(dat$male2N)),]
# remove empty rows
bsz <- bsz[!is.na(bsz$wingspan.mm.),]
bsz <- bsz[!is.na(bsz$hostType),]
# remove duplicated species
bsz.uniqe.sp <- unique(bsz$speciesNames)
b.size <- matrix(data = NA,
                 nrow = length(bsz.uniqe.sp),
                 ncol = 3)
colnames(b.size) <- c("species", "wingspan", "hostType")
b.size <- as.data.frame(b.size)
for(i in 1:nrow(b.size)){
  b.size$species[i] <- bsz.uniqe.sp[i]
  b.size$wingspan[i] <- mean(bsz$wingspan.mm.[bsz$speciesNames == bsz.uniqe.sp[i]])
  b.size$hostType[i] <- sample(bsz$hostType[bsz$speciesNames == bsz.uniqe.sp[i]], 1)
}
# get the wingspam and host data to the chroms dataset
dat$wingspan <- NA
dat$host <- NA
for(i in 1:nrow(dat)){
  if(length(b.size$wingspan[b.size$species == dat$SpecisName[i]]) != 0){
    dat$wingspan[i] <- b.size$wingspan[b.size$species == dat$SpecisName[i]]
    dat$host[i] <- b.size$hostType[b.size$species == dat$SpecisName[i]]
  }
}
dat <- dat[!is.na(dat$wingspan),]
# get the overlap between the tree and trait data
overlap <- sp.matches(dat, trees)
dat.new <- overlap$chroms
tree.names <- overlap$name_corrections
# rename tree and hosts species names
tree <- trees[[sample(1:10,1)]]
tree <- keep.tip(tree, overlap$name_corrections$name_on_tree)
# rename tree tip and host species name to reflect genera level matches
for(i in 1:Ntip(tree)){
  tree$tip.label[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == tree$tip.label[i]]
}
# sort data to match with the phylogenetic order
dat.sorted <- as.data.frame(matrix(data =NA,
                                   ncol = ncol(dat.new),
                                   nrow = nrow(dat.new)))
colnames(dat.sorted) <- colnames(dat.new)
for(i in 1:nrow(dat.new)){
  dat.sorted[i,] <- dat.new[dat.new$SpecisName == tree$tip.label[i],]
}
phy <- tree
phy.dat <- dat.sorted
rownames(phy.dat) <- phy.dat$SpecisName
# remove unwanted objects from the environment
rm(list = ls()[-c(10,11)])
phylo.lmod <- phylolm(male2N ~ wingspan : host, 
                     data = phy.dat, 
                     phy = phy,
                     model = "BM")
summary(phylo.lmod)
rm(list = ls())