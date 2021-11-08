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

# read poor runs
poorMCMC <- read.csv("../results/5a.fase.postivie.poorMCMC.csv")
poorMCMC <- poorMCMC$x
poorMCMC <- sort(poorMCMC)


# pars

# define number of clusters for parallel computing
NumberOfClusters <- length(poorMCMC)
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)

# define parameters for MCMC
iter.temp <- 20
tree.rep <- rep(1:10, each = 10)
iter <- 100
upper <- 200
prior <- make.prior.exponential(r = 5)
results <- vector(mode = "list", length = 100)
tree.depth <- vector(mode = "numeric", length = 100)



# read data
dat <- read.csv("../results/2.Lepidoptera.rate.analysis.feeding.csv")

x <- foreach(j = poorMCMC, .verbose = T, .packages = c("ape","maps","diversitree", "phytools", "chromePlus")) %dopar% {
  
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
  # get the mean transition rate
  trans12 <- mean(dat$tran12)
  trans21 <- mean(dat$tran21)
  # make the model
  model <- matrix(data = 0,
                  ncol = 2,
                  nrow = 2)
  colnames(model) <- rownames(model) <- c(0,1)
  # fill the model
  model[2,1] <- trans12
  model[2,2] <- -trans12
  model[1,1] <- -trans21
  model[1,2] <- trans21
  # get the ancestral state of the binary character
  root <- fitMk(tree = tree,
                x = setNames(chrom.sorted$feeding, chrom.sorted$binomial),
                fixedQ = model,
                pi = "fitzjohn")
  # simulate a binary trait
  sim.trait <-  sim.character(tree = tree,
                              pars = model,
                              model = "mkn",
                              x0 = as.numeric(sample(names(root$pi),
                                                     1, 
                                                     prob = root$pi))+1)
  # simulated character should be a binary character
  # perform a check to see if the simulated character is binary
  if(length(unique(sim.trait)) != 2){
    while (length(unique(sim.trait)) != 2) {
      sim.trait <-  sim.character(tree = tree,
                                  pars = model,
                                  model = "mkn",
                                  x0 = as.numeric(sample(names(root$pi),
                                                         1, 
                                                         prob = root$pi))+1)
    }
  }
  # format data to MCMC analysis
  # make tree unit length
  # make
  tree.depth[j] <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth[j]
  # make a new data table which has species names and chromosome number
  MCMC.dat <- as.data.frame(matrix(data = NA,
                                   nrow = nrow(chrom.sorted),
                                   ncol = 3))
  colnames(MCMC.dat) <- c("Species", "Chroms","simTrait")
  MCMC.dat$Species <- chrom.sorted$binomial
  MCMC.dat$Chroms <- chrom.sorted$haploid
  for(i in 1:nrow(MCMC.dat)){
    if(sim.trait[MCMC.dat$Species[i]] == 1){
      MCMC.dat$simTrait[i] <- 1
    }
    if(sim.trait[MCMC.dat$Species[i]] == 2){
      MCMC.dat$simTrait[i] <- 0
    }
  }
  # get the range of chromosome number
  rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
           range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
  # make a probability matrix for chromosome number
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
  con.lik <- constrainMuSSE.custom(data = chrom.mat,
                                   lik = lik,
                                   polyploidy = F,
                                   hyper = T,
                                   s.lambda = T,
                                   s.mu = T,
                                   oneway = F,
                                   constrain = list(drop.demi = T,
                                                    drop.poly = F))
  if(j %% 10 == 0){
    print(paste("processing replicate: ", 10 ," of tree ",tree.rep[j],
                sep = ""))
  }else{
    print(paste("processing replicate: ", (j %% 10) ," of tree ",tree.rep[j],
                sep = ""))
  }
  # make a temp MCMC to get the w parameter
  temp <- c()
  temp <- mcmc(lik = con.lik,
               x.init = runif(min=0, max=1,
                              n=length(argnames(con.lik))),
               prior = prior,
               nsteps = iter.temp,
               w = 1, 
               upper = upper)
  # get values for w
  w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
  # run MCMC
  results[[j]] <- mcmc(lik = con.lik,
                       x.init = runif(min=0, max=1,
                                      n=length(argnames(con.lik))),
                       nsteps = iter,
                       w = w,
                       prior = prior,
                       upper = upper)
  
  # convert rate parameter for millions of years
  results[[j]] <- results[[j]] / tree.depth[j]
}
# stop the cluster
stopCluster(cl)
# save
save.image("../results/5a.false.positive.test.Lepidoptera.rate.analysis.feeding.fix.incomplete.RData")
