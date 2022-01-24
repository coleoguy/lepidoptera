# Terrence Sylvester
# 28 August 2021
# pradakshanas@gmail.com

# load required libraries
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)
library(castor)

# load helper functions
source("helper.functions.R")

# define parameters for MCMC
iter.temp <- 20
tree.rep <- rep(1:10, each = 10)
iter <- 100
upper <- 200
prior <- make.prior.exponential(r = 5)
results <- vector(mode = "list", length = 100)

# load data
post.burnin <- read.csv("../results/2d.Lepidoptera.rate.analysis.proc.csv")

# read in data
trait <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
chrom <- read.delim("../data/chroms/chroms.txt", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# set J
j <- sample(1:length(trees),size = 1)

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
# get the overlap between the species and chromosome number data sets
sp.overlap <-  sp.matches(dat = chrom.dedup, trees = trees, hyper = T)
# processed dat
processed.dat <- sp.overlap$chroms
# rename tree tip names to reflect genera level matches
tree <- trees[[j]]
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
chrom.mat <- datatoMatrix(x = MCMC.dat[,c(1,2)],
                          range = rng,
                          hyper = F)
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
                        hyper = F,
                        # s.lambda = T,
                        # s.mu = T,
                        oneway = F,
                        verbose = T,
                        constrain = list(drop.demi = T,
                                         drop.poly = F))
# get the qmat mat
parMat <- con.lik[[2]]
# here there are rates that we are not interested at
# discard them from the table
parMat[!(parMat %in% c(1,2,5))] <- 0

# rate1 ascending aneuploidy - diploid or state 1 of hypertrait
# rate2 descending aneuploidy - diploid or state 1 of hypertrait
# rate5 polyploidization of a diploid or state 1 of hypertrait

# fill the qmat
parMat[parMat == 1] <- mean(post.burnin$asc1) * tree.depth
parMat[parMat == 2] <- mean(post.burnin$desc1) * tree.depth
parMat[parMat == 5] <- mean(post.burnin$pol1) * tree.depth

# rename the columns so that states now start from 1
states <- 1:ncol(chrom.mat)
names(states) <- colnames(chrom.mat)
colnames(parMat) <- rownames(parMat) <- colnames(chrom.mat) <- states

# make row sum to zero
for(i in 1:nrow(parMat)){
  parMat[i,i] <- 0 - sum(parMat[i,])
}

# sort chrom mat to match with the tree tip order
chrom.mat <- chrom.mat[tree$tip.label,]
# sort MCMC dat to match with the tree tip order
rownames(MCMC.dat) <- MCMC.dat$Species
MCMC.dat <- MCMC.dat[tree$tip.label,]

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

# get chromosome number at the tip
dat.tips <- as.data.frame(matrix(data = NA,
                                 nrow = Ntip(tree),
                                 ncol = 2))
colnames(dat.tips) <- c("tip", "chrom")
for(i in 1:Ntip(tree)){
  dat.tips$tip[i] <- i
  dat.tips$chrom <- MCMC.dat$Chroms[MCMC.dat$Species == tree$tip.label[i]]
}
# get the branch specific rates
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

# save the results
save.image("../results/9a.branch.specific.rates.Lepidoptera.RData")
