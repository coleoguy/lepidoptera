# Terrence Sylvester
# pradakshanas@gmail.com
# 24 Jan 2021

library(diversitree)
library(ape)
library(phytools)
library(chromePlus)

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
                          constrain = list(drop.demi = T,
                                           drop.poly = F),
                          s.lambda = T,
                          s.mu = T)
# run a mcmc for smaller number of generations to get the
# w parameter
temp <- mcmc(lik = con.lik,
             x.init = runif(min=0, max=1,
                            n=length(argnames(con.lik))),
             nsteps = iter.temp,
             w = 1,
             prior = prior,
             upper = upper)
# get values for w
w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik)) + 1)], quantile, c(.05, .95)))
# run the mcmc
results.emp <- mcmc(lik = con.lik,
                    x.init = runif(min=0, max=1,
                                   n=length(argnames(con.lik))),
                    nsteps = iter,
                    w = w,
                    prior = prior,
                    upper = upper)
# convert rate parameter for millions of years
# results.emp <- results.emp[2:(length(argnames(con.lik)) + 1)] / tree.depth

# discard the burnin portion
# lets discard these trees from our resutls
proc.results <- results.emp

# sample the burnin portion
# burnin 
burnin <- .50
burn <- -1:-(iter*burnin)
post.burnin <- results.emp[burn,]

# get likelihood after dropping a tip
tip.names <- tree$tip.label
# make a table to store likelihoods
lik.dat <- as.data.frame(matrix(data = NA,
                                nrow = length(tip.names),
                                ncol = 3))
colnames(lik.dat) <- c("tip.number", "tip.name", "likelihood")
# fill in rest of the datatable
lik.dat$tip.number[i] <- 1:length(tip.names)
lik.dat$tip.name[i] <- tip.names

# get the likelihood of the complete tree
tree.lik <- con.lik(colMeans(post.burnin[,2:(length(argnames(con.lik)) + 1)]))

# calculate the likelihood by dropping a tip
for(i in 1:length(tip.names)){
  # drop a single tip
  sub.tree <- drop.tip(tree, tip.names[i])
  # remove that data from the data table
  sub.MCMC.dat <- MCMC.dat[MCMC.dat$Species != tip.names[i],]
  # make the likelihood function
  # get the range of chromosome number
  sub.rng <- c(range(sub.MCMC.dat$Chroms, na.rm = T)[1] - 1,
               range(sub.MCMC.dat$Chroms, na.rm = T)[2] + 1)
  # our states are coded such that specialists will be state 1 and
  # generalists will be state 2
  sub.chrom.mat <- datatoMatrix(x = sub.MCMC.dat,
                                range = sub.rng,
                                hyper = F)
  # make the likelihood function
  sub.lik <- make.musse(tree = sub.tree,
                        states = sub.chrom.mat,
                        k = ncol(sub.chrom.mat),
                        strict = FALSE,
                        control = list(method = "ode"))
  # constrain the likelihood function
  sub.con.lik <- constrainMuSSE(data = sub.chrom.mat,
                                lik = sub.lik,
                                polyploidy = F,
                                hyper = F,
                                constrain = list(drop.demi = T,
                                                 drop.poly = F),
                                s.lambda = T,
                                s.mu = T)
  # get likelihood
  lik.dat$likelihood[i] <-  sub.con.lik(colMeans(post.burnin[,2:(length(argnames(con.lik)) + 1)]))
}
save.image("../results/drop.taxa.and.calculate.likelihood.RData")
