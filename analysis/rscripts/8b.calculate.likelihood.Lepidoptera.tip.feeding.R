# Terrence Sylvester
# 28 August 2021
# pradakshanas@gmail.com

# load required libraries
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)
library(doSNOW)

# load helper functions
source("helper.functions.R")

dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv", as.is = T)

# define number of clusters for parallel computing
NumberOfClusters <- 4
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)

# define parameters for MCMC
iter.temp <- 20
tree.rep <- rep(1:10, each = 10)
iter <- 100
upper <- 200
prior <- make.prior.exponential(r = 5)
results <- vector(mode = "list", length = 100)
# tree.depth <- vector(mode = "numeric", length = 100)

# set J
j <- 1

# read in data
trait <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
chrom <- read.delim("../data/chroms/chroms.txt", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# fill in binomial column so that it does not includes sub species names
for(i in 1:nrow(chrom)){
  chrom$binomial[i] <- paste(chrom$Genus[i], chrom$species[i], sep = "_")
}

# filter out unwanted columns
chrom <- chrom[c("Family","Genus","binomial","female2N", "male2N", "notesMale2N")]

# process data
# sample a single chromosome number when there is a range
proc.chrom <- chromSampler(dat = chrom)
# remove species that lack male chromosome number data
proc.chrom <-  proc.chrom[!(is.na(proc.chrom$male2N)),]
# remove unwanted columns from the data table
proc.chrom <- proc.chrom[c("Family","Genus","binomial","male2N")]
# remove species that do not have information on chromosome number
proc.chrom <- proc.chrom[!is.na(proc.chrom$male2N),]
# get the haploid chromosome number
proc.chrom$haploid <- NA
proc.chrom$haploid <- proc.chrom$male2N / 2

# if there are species with floating point values for haploid number. 
# randomly sample a single value (ceiling or floor) for these species
for(i in 1:nrow(proc.chrom)){
  if(proc.chrom$haploid[i] %% 1 != 0){
    proc.chrom$haploid[i] <- sample(c(ceiling(proc.chrom$haploid[i]), floor(proc.chrom$haploid[i])),1)
  }
}
# check if there are species with multiple records for chromosome number
# if present, randomly sample a single value
chrom.dedup <- as.data.frame(matrix(data = NA,
                                    nrow = length(unique(proc.chrom$binomial)),
                                    ncol = ncol(proc.chrom)))
colnames(chrom.dedup) <- colnames(proc.chrom)

# fill out species names
chrom.dedup$binomial <- unique(proc.chrom$binomial)
# use a for loop to fill out rest of the information
# for species that have multiple records we first store them in a temporary object
# and then sample a single value
for (i in 1:nrow(chrom.dedup)) {
  if(sum(proc.chrom$binomial == chrom.dedup$binomial[i]) == 1){
    chrom.dedup$Family[i] <- proc.chrom$Family[proc.chrom$binomial == chrom.dedup$binomial[i]]
    chrom.dedup$Genus[i] <- proc.chrom$Genus[proc.chrom$binomial == chrom.dedup$binomial[i]]
    chrom.dedup$haploid[i] <- proc.chrom$haploid[proc.chrom$binomial == chrom.dedup$binomial[i]]
    
  }
  if(sum(proc.chrom$binomial == chrom.dedup$binomial[i]) > 1){
    hit <- proc.chrom$haploid[proc.chrom$binomial == chrom.dedup$binomial[i]] 
    chrom.dedup$haploid[i] <- sample(hit, 1)
  }
  chrom.dedup$Family[i] <- sample(proc.chrom$Family[proc.chrom$binomial == chrom.dedup$binomial[i]],1)
  chrom.dedup$Genus[i] <- sample(proc.chrom$Genus[proc.chrom$binomial == chrom.dedup$binomial[i]],1)
}
# remove unwanted columns
chrom.dedup <- chrom.dedup[,-4]
# get the development state of each species
chrom.dedup$feeding <- NA
for(i in 1:nrow(chrom.dedup)){
  hit.feeding <- trait$type[trait$binomial == chrom.dedup$binomial[i]]
  if(length(hit.feeding) != 0){
    chrom.dedup$feeding[i] <- trait$type[trait$binomial == chrom.dedup$binomial[i]]
  }
}
# get the genera level representations for development
# get genera which have no development data
gen.feeding.no <- chrom.dedup$Genus[is.na(chrom.dedup$feeding)]
# get genera which have development data
gen.feeding.yes <- chrom.dedup$Genus[!is.na(chrom.dedup$feeding)]
# get genera that are absent the subset of genera that have development data
gen.feeding.no.match <- c()
for(i in 1:length(gen.feeding.no)){
  if(gen.feeding.no[i] %in% gen.feeding.yes == F){
    gen.feeding.no.match <- c(gen.feeding.no.match, gen.feeding.no[i])
  }
}
# get the unique list of names of these genera
gen.feeding.no.match.unique <- unique(gen.feeding.no.match)
# fill these genera with development data
for(i in 1:length(gen.feeding.no.match.unique)){
  if(gen.feeding.no.match.unique[i] %in% chrom.dedup$Genus){
    hit.feeding.genera <- trait$type[trait$genus == gen.feeding.no.match.unique[i]]
    if(length(hit.feeding.genera) != 0){
      if(length(unique(hit.feeding.genera)) == 1){
        chrom.dedup$feeding[chrom.dedup$Genus == gen.feeding.no.match.unique[i]] <- unique(hit.feeding.genera)
      }else{
        chrom.dedup$feeding[chrom.dedup$Genus == gen.feeding.no.match.unique[i]] <- sample(hit.feeding.genera,1)
      }
    }
  }
}
# remove species that do not have Development information
chrom.dedup <- chrom.dedup[!is.na(chrom.dedup$feeding),]
# get the overlap between the species and chromosome number data sets
sp.overlap <-  sp.matches(dat = chrom.dedup, trees = trees, hyper = T)
# processed dat
processed.dat <- sp.overlap$chroms
# rename tree tip names to reflect genera level matches
tree <- trees[[tree.rep[j]]]
tree <- keep.tip(tree, sp.overlap$name_corrections$name_on_tree)
for(i in 1:Ntip(tree)){
  tree$tip.label[i] <- sp.overlap$name_corrections$species_name[sp.overlap$name_corrections$name_on_tree == tree$tip.label[i]]
}
# sort the chromosome number dataset so that it matches with the tip order
chrom.sorted <- as.data.frame(matrix(data = NA,
                                     nrow = nrow(processed.dat),
                                     ncol = ncol(processed.dat)))
colnames(chrom.sorted) <- colnames(processed.dat)
for(i in 1:nrow(chrom.sorted)){
  chrom.sorted[i,] <- processed.dat[which(processed.dat$binomial == tree$tip.label[i]),]
}
# format data to MCMC analysis
# make tree unit length
# make
tree.depth <- max(branching.times(tree))
tree$edge.length <- tree$edge.length / tree.depth
# make a new data table which has species names and chromosome number
MCMC.dat <- as.data.frame(matrix(data = NA,
                                 nrow = nrow(chrom.sorted),
                                 ncol = 3))
colnames(MCMC.dat) <- c("Species", "Chroms","feeding")
MCMC.dat$Species <- chrom.sorted$binomial
MCMC.dat$Chroms <- chrom.sorted$haploid
MCMC.dat$feeding <- chrom.sorted$feeding
# in the dataset of larval feeding 1 means specialists and 0 means generalists
# change that so 1 means generalist and 0 means specialist
MCMC.dat$gen.prob <- NA
for(i in 1:nrow(MCMC.dat)){
  if(MCMC.dat$feeding[i] == 1){
    MCMC.dat$gen.prob[i] <- 0
  }
  if(MCMC.dat$feeding[i] == 0){
    MCMC.dat$gen.prob[i] <- 1
  }
}
# get the range of chromosome number
rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
         range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
# make a probability matrix for chromosome number
chrom.mat <- datatoMatrix(x = MCMC.dat[,c(1,2,4)],
                          range = rng,
                          hyper = T)
# make the likelihood function
lik <- make.mkn(tree = tree,
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = FALSE,
                control = list(method = "ode"))
# constrain the likelihood function
con.lik <- constrainMkn(data = chrom.mat,
                        lik = lik,
                        polyploidy = F,
                        hyper = T,
                        # s.lambda = T,
                        # s.mu = T,
                        oneway = F,
                        constrain = list(drop.demi = T,
                                         drop.poly = F))

# get likelihood after dropping a tip
tip.names <- tree$tip.label
# make a table to store likelihoods
lik.table <- as.data.frame(matrix(data = NA,
                                  nrow = length(tip.names),
                                  ncol = 3))
colnames(lik.table) <- c("tip.number", "tip.name", "likelihood")
# fill in rest of the datatable
lik.table$tip.number <- 1:length(tip.names)
lik.table$tip.name <- tip.names
# get the likelihood of the complete tree
full.tree.lik <- con.lik(colMeans(dat[,2:(length(argnames(con.lik)) + 1)])*tree.depth)

# calculate the likelihood by dropping a tip
x <- foreach(i = length(tree$tip.label), .verbose = T, .packages = c("ape","maps","diversitree", "phytools", "chromePlus")) %dopar% {
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
  sub.chrom.mat <- datatoMatrix(x = sub.MCMC.dat[,c(1,2,4)],
                                range = sub.rng,
                                hyper = T)
  # make the likelihood function
  sub.lik <- make.mkn(tree = sub.tree,
                        states = sub.chrom.mat,
                        k = ncol(sub.chrom.mat),
                        strict = FALSE,
                        control = list(method = "ode"))
  # constrain the likelihood function
  sub.con.lik <- constrainMkn(data = sub.chrom.mat,
                              lik = sub.lik,
                              polyploidy = F,
                              hyper = T,
                              # s.lambda = T,
                              # s.mu = T,
                              oneway = F,
                              constrain = list(drop.demi = T,
                                               drop.poly = F))
  # get likelihood
  lik.table$likelihood[i] <-  sub.con.lik(colMeans(dat[,2:(length(argnames(sub.con.lik)) + 1)])*tree.depth)
}
# stop the cluster
stopCluster(cl)
# save
save.image("../results/8b.calculate.likelihood.Lepidoptera.tip.feeding.RData")
