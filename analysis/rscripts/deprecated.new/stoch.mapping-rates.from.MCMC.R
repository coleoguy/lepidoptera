# Terrence Sylvester
# 23 February 2021
# pradakshanas@g,ail.com

# load libraries
library(phytools)
library(diversitree)

# get helper functions
source("../rscripts/helper.functions.R")

# read in data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# combine genus and species names to get the full name in hosts dataset 
host$binomial <- paste(host$genus, host$species, sep = "_")

# read in post burnin portion of the mcmc 
post.burnin <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
# mean genaralist to specialist rate
q01 <- mean(post.burnin$tran21) # in mcmc state 1 was specialist and 2 was generalist
# mean specialist to generalist rate
q10 <- mean(post.burnin$tran12)

# place holder for simmap results
simmap.MCMC.rt <- simmap.MCMC <- vector(mode = "list", length = 10)
names(simmap.MCMC) <- names(simmap.MCMC.rt) <- paste("tree", 1:10, sep = ".")

# sample a single chromosome number when there is a range
new.dat <- chromSampler(dat = dat)
# remove species that lack male chromosome number data
new.dat <-  new.dat[!(is.na(new.dat$male2N)),]
# get the overlap between the traits and the phylogeny
overlap <- sp.matches(dat = new.dat, trees = trees, use.sub.species = F)
# processed dat
processed.dat <- overlap$chroms
# rename tree and hosts species names
phy <- vector(mode = "list", length = 10)
for(i in 1:10){
  phy[[i]] <- keep.tip(trees[[i]], overlap$name_corrections$name_on_tree)
}
class(phy) <- "multiPhylo"

qmat <-  matrix(data = c(-q01, q01, q10, -q10),
                nrow = 2,
                ncol = 2,
                byrow = T,dimnames = list(c(0,1),c(0,1)))

# rename tree tip and host species name to reflect genera level matches
for(j in 1:10){
  for(i in 1:Ntip(phy[[j]])){
    phy[[j]]$tip.label[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == phy[[j]]$tip.label[i]]
  }
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

mcmc.tip.states <- MCMC.dat$HostType
names(mcmc.tip.states) <- MCMC.dat$Species

# get simmaps on a random trait
for(i in 1:10){
  root.state <- sample(c(0,1),1)
  # get the tip states of a neutrally evolving trait
  tip.states <- sim.character(tree = phy[[i]],
                              pars = c(q01,q10), #(0 to 1, 1 to 0)
                              model = "mk2", 
                              x0 = root.state) # root state
  if(root.state == 0){
    pi <- c(1,0)
  }
  if(root.state == 1){
    pi <- c(0,1)
  }
  
  simmap.MCMC.rt[[i]] <- make.simmap(tree = phy[[i]],
                                  x = tip.states,
                                  Q = qmat,
                                  pi = pi,
                                  nsim = 100)
}

# get simmaps for the actual trait
simmap.MCMC <- make.simmap(tree = phy,
                           x = mcmc.tip.states,
                           Q = qmat,
                           pi = c(1,0),
                           nsim = 100)
# save
save.image("../results/stoch.mapping-rates-from-MCMC.RData")
