# load libraries
library(castor)
library(ape)
library(chromePlus)
library(viridis)
library(phytools)
library(diversitree)

# get helper functions
source("helper.functions.R")
source("tree.painter.R")

# read in data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# read in rates
post.burn.in <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)

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
# rename tree and hosts species names
tree <- trees[[sample(1:10,1)]]
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
# tree is not utrametric (eventhough it is a BEAST output). Make that correction
tree <- force.ultrametric(tree)
# remove unwanted tips
tree <- drop.tip(tree, processed.dat$binomial[is.na(processed.dat$male2N)])
# make tree unit length
tree.depth <- max(branching.times(tree))
tree$edge.length <- tree$edge.length / tree.depth
# remove unwanted hosts data
host <- host[host$binomial %in% tree$tip.label,]
# get the range of chromosome number
rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
         range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
# our states are coded such that specialists will be state 1 and
# generalists will be state 2
chrom.mat <- datatoMatrix(x = MCMC.dat,
                          range = rng,
                          hyper = T)
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
                          hyper = T,
                          constrain = list(drop.demi = T,
                                           drop.poly = F),
                          s.lambda = F,
                          s.mu = F,
                          verbose = T)
rate.mat <- con.lik[[2]]
rate.mat[!(rate.mat %in% c(1,2,3,4,5,6,8,9))] <- 0

rate.mat[rate.mat == 1] <- mean(post.burn.in$asc1)
rate.mat[rate.mat == 2] <- mean(post.burn.in$desc1)
rate.mat[rate.mat == 3] <- mean(post.burn.in$asc2)
rate.mat[rate.mat == 4] <- mean(post.burn.in$desc2)
rate.mat[rate.mat == 5] <- mean(post.burn.in$pol1)
rate.mat[rate.mat == 6] <- mean(post.burn.in$pol2)
rate.mat[rate.mat == 8] <- mean(post.burn.in$tran12)
rate.mat[rate.mat == 9] <- mean(post.burn.in$tran21)

for(i in 1:ncol(rate.mat)){
  rate.mat[i,i] <- -sum(rate.mat[i,])
}

# set the rate class of all the branches class with the multiplication 
# facter of 1

# set the rates and rate classes
rate.class <- c(1,2,3,4,5,6,7,8,9,10,11)
rate <- c(.1,.125,.1666,25,.5,1,2,4,6,8,10)

tree$rates <- rep(which(rate == 1),Nedge(tree))

# get the states of the tree
states <- 1:ncol(chrom.mat)
names(states) <- colnames(chrom.mat)

tip.states <- vector(mode = "numeric", length = Ntip(tree))

for(i in 1:Ntip(tree)){
  hit.row <- which(rownames(chrom.mat) == tree$tip.label[i])
  hit.col <- which(chrom.mat[hit.row,] == 1)
  tip.states[i] <-  as.numeric(states[hit.col])
}

colnames(rate.mat) <- rownames(rate.mat) <- states

painted.tree <- tree.paintR.ver.2(tree = tree,
                            tip_states = tip.states,
                            qmat = rate.mat,
                            iter = 1,
                            rate.class = rate.class,
                            rate = rate,
                            sample_mode = "median",
                            return_br_table = T)

plot(tree,
     show.tip.label = F,
     edge.color = viridis(n = length(rate.class))[painted.tree$tree$rates],
     edge.width = 2,type = "fan")




