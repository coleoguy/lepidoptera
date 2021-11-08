# Terrence Sylvester
# 28 August 2021
# pradakshanas@gmail.com

# load required libraries
library(ape)
library(chromePlus)
library(diversitree)
library(phytools)

# load helper functions
source("helper.functions.R")

# load data
dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.proc.csv")
# pars
tree.rep <- rep(1:10, each = 10)
# in mcmc state 1 was non direct development and 2 was direct development
# mean non direct developing to direct developing rate
q01 <- mean(dat$tran12) 
# mean direct developing to non direct developing rate
q10 <- mean(dat$tran21)
# mean rate of chromosome fission at state 1 
state.01.asc <- mean(dat$asc1)
# mean rate of chromosom fission at state 2
state.02.asc <- mean(dat$asc2)
# mean rate of chromosome fusion at state 1
state.01.desc <- mean(dat$desc1)
# mean rate of chromosome fusuin at state 2
state.02.desc <- mean(dat$desc2)
# mean rate of polyploidy at state 1
state.01.pol <- mean(dat$pol1)
# mean rate of polyploidy at state 2
state.02.pol <- mean(dat$pol2)

# place holder for final results
results.final <- vector(mode = "list", length = 100)

for(j in 1:100){
  # place holder for results
  results.chroms <- vector(mode = "list", length = 100)
  results.states <- vector(mode = "list", length = 100)
  results.temp <- c()
  
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
  # make states
  chrom.states <- as.data.frame(matrix(data = NA,
                                       nrow = length(rng[1]:rng[2]) * 2,
                                       ncol = 3))
  colnames(chrom.states) <- c("chroms", "binary", "state")
  # fill chrom states
  chrom.states$chroms <- rep(rng[1]:rng[2], 2)
  chrom.states$binary <- rep(c(1,0), each = length(rng[1]:rng[2]))
  
  for(i in 1:nrow(chrom.states)){
    if (nrow(chrom.states) < 100) 
      pad <- 2
    if (nrow(chrom.states) >= 100) 
      pad <- 3
    if (nrow(chrom.states) < 10) 
      pad <- 1
    chrom.states$state[i] <- sprintf(paste("%0", pad, "d", sep = ""),i)
  }
  
  MCMC.dat$State <- NA
  # properly name each state in the qmatrix
  for(i in 1:nrow(MCMC.dat)){
    hit.chrom <- chrom.states$chroms == MCMC.dat$Chroms[i]
    hit.binary <- chrom.states$binary == MCMC.dat$feeding[i]
    hit.state<- which(hit.binary & hit.chrom)
    MCMC.dat$State[i] <- chrom.states$state[hit.state]
  }
  # build qmatrix
  # get pars
  pars = c(state.01.asc,   # ascending at state 1
           state.02.asc,   # ascending at state 2
           state.01.desc,  # descending at state 1
           state.02.desc,  # descending at state 2
           0,              # demiploidy at state 1
           0,              # demiploidy at state 2
           state.01.pol,   # polyploidy at state 1
           state.02.pol,   # polyploidy at state 2
           q01,            # transition from state 1 to state 2
           q10)            # transition from state 2 to state 1
  q <- buildQ(pars = pars, limits = rng)
  # get the root of the tree
  for(i in 1:nrow(q)){
    q[i,i] <- -sum(q[i,])
  }
  # rename the colnames and the row names of the qmatrix
  colnames(q) <- rownames(q) <- chrom.states$state
  # make two copies of the qmatrix to be used in two different functions
  q.model.fit <- q
  q.sim.char <- q  
  # prepare inputs for fitMK
  colnames(q.model.fit) <- rownames(q.model.fit) <- colnames(chrom.mat)
  # get root states
  model.fit <- fitMk(tree = tree,
                     x = chrom.mat,
                     fixedQ = q.model.fit, 
                     pi = "fitzjohn")
  # rename the col names of q matrix used in simulated datasets
  colnames(q.sim.char) <- rownames(q.sim.char) <- 1:ncol(q.sim.char)
  names(model.fit$pi) <- 1:ncol(q.sim.char)
  # simulate datasets
  for(k in 1:100){
    print(paste("Tree", j, "simulation", k))
    # simulate datasets
    sim.chroms <- sim.character(tree = tree,
                                pars = q.sim.char,
                                x0 = as.numeric(sample(names(model.fit$pi),size = 1,prob = model.fit$pi)),
                                model = "mkn")
    # get states
    states <- vector(mode = "numeric", length = length(sim.chroms))
    midpoint <-  chrom.states$state[which(chrom.states$binary == 1 & chrom.states$chroms == max(chrom.states$chroms))]
    states[sim.chroms > as.numeric(midpoint)] <- 0
    states[sim.chroms <= as.numeric(midpoint)] <- 1
    # get chroms
    sim.chrom.states <- vector(mode = "numeric", length = length(sim.chroms))
    for(i in 1:length(sim.chroms)){
      sim.chrom.states[i] <-  chrom.states$chroms[as.numeric(chrom.states$state) == sim.chroms[i]]
    }
    results.chroms[[k]] <- sim.chrom.states
    results.states[[k]] <- states
  }
  results.temp <- list(results.chroms,
                       results.states)
  names(results.temp) <- c("chroms",
                           "states")
  results.final[[j]] <- results.temp
}

save.image("../results/6a.model.adequacy.test.Lepidoptera.rate.analysis.feeding.RData")
