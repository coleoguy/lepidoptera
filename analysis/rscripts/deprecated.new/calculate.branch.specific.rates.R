# Terrence Sylvester
# pradakshanas@gmail.com
# 24 Jan 2021

library(diversitree)
library(ape)
library(phytools)
library(chromePlus)
library(castor)
library(viridis)

# read in rates
load("../results/drop.taxa.and.calculate.likelihood.RData")

# remove all unwanted data
rm(list = ls()[-19])

# get helper functions
source("../rscripts/helper.functions.R")

# parameters for MCMC
iter <- 100                                       
iter.temp <- 20                 
tree.num <- rep(1:10, each=10)                    
prior <- make.prior.exponential(r=2)              
upper <- 50
results <- vector(mode = "list", length = 100)

# run the MCMC
# read in data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")
# combine genus and species names to get the full name in hosts dataset and in
# chroms dataset
host$binomial <- paste(host$genus, host$species, sep = "_")
# sample a single chromosome number when there is a range
new.dat <- chromSampler(dat = dat)
# remove species that lack male chromosome number data
new.dat <-  new.dat[!(is.na(new.dat$male2N)),]
# get the overlap between the traits and the phylogeny
overlap <- sp.matches(dat = new.dat, trees = trees, use.sub.species = F)
# processed dat
processed.dat <- overlap$chroms
# pick a single tree
tree <- trees[[1]]
tree <- keep.tip(tree, overlap$name_corrections$name_on_tree)
# rename tree tip and host species name to reflect genera level matches
for(i in 1:Ntip(tree)){
  tree$tip.label[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == tree$tip.label[i]]
}
host.new <- host[host$binomial %in% overlap$name_corrections$name_on_tree,]
for(i in 1:nrow(host.new)){
  host.new$binomial[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == host.new$binomial[i]]
}
# make a new data table which has species names chromosome number and host type
MCMC.dat <- as.data.frame(matrix(data = NA,
                                 nrow = nrow(processed.dat),
                                 ncol = 3))
colnames(MCMC.dat) <- c("Species", "Chroms", "HostType")
MCMC.dat$Species <- processed.dat$binomial
MCMC.dat$Chroms <- processed.dat$male2N / 2
for(i in 1:nrow(MCMC.dat)){
  MCMC.dat$HostType[i] <- host.new$type[host.new$binomial == MCMC.dat$Species[i]]
}
# tree is not utrametric (even though it is a BEAST output). Make that correction
tree <- force.ultrametric(tree)
# remove unwanted tips
tree <- drop.tip(tree, processed.dat$binomial[is.na(processed.dat$male2N)])
# make tree unit length
tree.depth <- max(branching.times(tree))
tree$edge.length <- tree$edge.length / tree.depth
# remove unwanted hosts data
host <- host[host$SpeciesNames %in% tree$tip.label,]
# get the range of chromosome number
rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
         range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
# our states are coded such that specialists will be state 1 and
# generalists will be state 2
chrom.mat <- datatoMatrix(x = MCMC.dat,
                          range = rng,
                          hyper = F)
# make the likelihood function
lik <- make.musse(tree = tree,
                  states = chrom.mat,
                  k = ncol(chrom.mat),
                  strict = FALSE,
                  control = list(method = "ode"))
# constrain the likelihood function
con.lik <- constrainMuSSE(data = chrom.mat,
                          lik = lik,
                          polyploidy = F,
                          hyper = F,
                          verbose = T,
                          constrain = list(drop.demi = T,
                                           drop.poly = F))
# get the par mat
parMat <- con.lik[[2]]
# here there are rates that we are not interested at
# discard them from the table
parMat[!(parMat %in% c(1,2,5))] <- 0

# rate1 ascending aneuploidy - diploid or state 1 of hypertrait
# rate2 descending aneuploidy - diploid or state 1 of hypertrait
# rate5 polyploidization of a diploid or state 1 of hypertrait

# fill the qmat
parMat[parMat == 1] <- mean(post.burnin$asc1)
parMat[parMat == 2] <- mean(post.burnin$desc1)
parMat[parMat == 5] <- mean(post.burnin$pol1)
# rename the columns so that states now start from 1
states <- 1:ncol(chrom.mat)
names(states) <- colnames(chrom.mat)
# sort chrom mat to match with the tree tip order
chrom.mat <- chrom.mat[tree$tip.label,]
# sort MCMC dat to match with the tree tip order
rownames(MCMC.dat) <- MCMC.dat$Species
MCMC.dat <- MCMC.dat[tree$tip.label,]
# make row sum to zero
for(i in 1:nrow(parMat)){
  parMat[i,i] <- 0 - sum(parMat[i,])
}
colnames(parMat) <- rownames(parMat) <- colnames(chrom.mat) <- states
# get the ancestral states
asr <-  asr_mk_model(tree = tree,
                     tip_states = NULL,
                     tip_priors = chrom.mat,
                     transition_matrix = parMat,
                     Nstates = ncol(parMat),
                     include_ancestral_likelihoods = T,
                     root_prior = "empirical")
# get the most likely chromosome number state for each node
node.chroms <- as.data.frame(matrix(data = NA,
                                    nrow = Nnode(tree),
                                    ncol = 3))
colnames(node.chroms) <- c("node", "state", "chromNum")
# fill the table
for(i in 1:nrow(node.chroms)){
  node.chroms$node[i] <- i + Ntip(tree)
  poss.state <- which.max(asr$ancestral_likelihoods[i,])
  if(length(poss.state) > 1){
    node.chroms$state[i] <- sample(poss.state,1) 
  }else{
    node.chroms$state[i] <- poss.state     
  }
  node.chroms$chromNum[i] <- as.numeric(names(states[states == node.chroms$state[i]]))
}

# get tip values
dat.tips <- as.data.frame(matrix(data = NA,
                                 nrow = Ntip(tree),
                                 ncol = 2))

colnames(dat.tips) <- c("tip", "chrom")

for(i in 1:Ntip(tree)){
  dat.tips$tip[i] <- i
  dat.tips$chrom <- MCMC.dat$Chroms[MCMC.dat$Species == tree$tip.label[i]]
}

branch.rates <- as.data.frame(matrix(data = NA,
                                     nrow = Nedge(tree),
                                     ncol = 3))



colnames(branch.rates) <- c('branch', "rate", "rateClass")
counter <- 1
for(i in 1:nrow(node.chroms)){
  hit.edge <- which(tree$edge[,1] == node.chroms$node[1])
  for(j in 1:length(hit.edge)){
    end.node <- tree$edge[hit.edge[j],2]
    if(end.node %in% c(1:Ntip(tree))){
      br.rate <- abs(dat.tips$chrom[end.node] - node.chroms$chromNum[i]) / tree$edge.length[hit.edge[j]]
      branch.rates$branch[counter] <- counter
      branch.rates$rate[counter] <- br.rate
      counter <- counter + 1
    }
    else{
      br.rate <- abs(node.chroms$chromNum[node.chroms$node == end.node] - node.chroms$chromNum[i]) / tree$edge.length[hit.edge[j]]
      branch.rates$branch[counter] <- counter
      branch.rates$rate[counter] <- br.rate
      counter <- counter + 1
    }
  }
}

# scale branch rates to millions of years
# branch.rates$rate <- branch.rates$rate / tree.depth

for(i in 1:nrow(branch.rates)){
  if(branch.rates$rate[i] <= quantile(x = branch.rates$rate, 0.25)){
    branch.rates$rateClass[i] <- 1
  }
  if(branch.rates$rate[i] > quantile(x = branch.rates$rate, 0.25) & branch.rates$rate[i] <= quantile(x = branch.rates$rate, 0.75)){
    branch.rates$rateClass[i] <- 2
  }
  if(branch.rates$rate[i] > quantile(x = branch.rates$rate, 0.75)){
    branch.rates$rateClass[i] <- 3
  }
}
# 
# plot(tree, type = "f", show.tip.label = F,
#      edge.color = viridis(3)[branch.rates$col], edge.width = 2)
# 
# text(x = 0.75,
#      y = 1.05+ 0.25,
#      labels = "Rate class",
#      pos = 4,
#      cex = 0.9)
# 
# points(x = rep(.8,3),
#        y = c(1.0,0.95,0.90),
#        pch = 16,
#        col = viridis(3))
# 
# text(x = rep(.8,3),
#      y = c(1.0,0.95,0.90),
#      labels = c("Low", "Medium", "High"),
#      pos = 4, cex = .8)
# 
# tiplabels(pch = 16, col = c("blue","red")[MCMC.dat$HostType + 1], cex = .5, offset = .02)
# 
# points(x = rep(.8,2),
#        y = c(0.79,0.74),
#        pch = 16,
#        col = c("blue","red"))
# 
# text(x = 0.75,
#      y = 0.84,
#      labels = "Binary trait",
#      pos = 4, cex = .9)
# 
# text(x = rep(.8,2),
#      y = c(0.79,0.74) + 0.25,
#      labels = c("Generalists", "Specialists"),
#      pos = 4,
#      cex = .8)

# make a maps object and add it to tree in order to colour branches
maps <- vector(mode = "list", length = Nedge(tree))
for(i in 1:length(maps)){
  maps[[i]] <- tree$edge.length[i]
  names(maps[[i]]) <- branch.rates$rateClass[i]
}
tree$maps <- maps

plotTree.wBars(tree = tree,
               x = setNames(MCMC.dat$Chroms, MCMC.dat$Species),
               type = "fan",
               col = c("blue","red")[MCMC.dat$HostType + 1],
               border = NA,
               color = setNames(c(viridis(3)),
                                c(1,2,3)),
               method = "plotSimmap")


text(x = 0.75,
     y = 1.05 + 0.25,
     labels = "Rate class",
     pos = 4,
     cex = 0.8)

points(x = rep(.8,3),
       y = c(1.0,0.95,0.90)+ 0.25,
       pch = 16,
       col = viridis(3))

text(x = rep(.8,3),
     y = c(1.0,0.95,0.90)+ 0.25,
     labels = c("Low", "Medium", "High"),
     pos = 4, cex = .7)

points(x = rep(.8,2),
       y = c(0.79,0.74)+ 0.23,
       pch = 16,
       col = c("blue","red"))

text(x = 0.75,
     y = 0.84+ 0.23,
     labels = "Binary trait",
     pos = 4, cex = .8)

text(x = rep(.8,2),
     y = c(0.79,0.74) + 0.23,
     labels = c("Generalists", "Specialists"),
     pos = 4,
     cex = 0.7)
