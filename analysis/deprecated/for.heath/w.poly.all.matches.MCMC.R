# Terrence Sylvester
# pradakshanas@gmail.com
# 24 Jan 2021

# load libraries
library(doSNOW)
# library(parallel)
library(diversitree)
library(ape)
library(phytools)
library(chromePlus)

# get helper functions
source("../rscripts/helper.functions.R")

# define number of clusters for parallel computing
NumberOfClusters <- 2
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)
# parameters for MCMC
iter <- 100                                       
iter.temp <- 20                 
tree.num <- rep(1:10, each=10)                    
prior <- make.prior.exponential(r=2)              
upper <- 50
results <- vector(mode = "list", length = 100)

# run the MCMC
x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","maps","diversitree", "phytools", "chromePlus")) %dopar% {
  print(1)
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
  # rename tree and hosts species names
  tree <- trees[[tree.num[j]]]
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
                          s.mu = F)
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
  w <- diff(sapply(temp[11:20, 2:9], quantile, c(.05, .95)))
  # run the mcmc
  results[[j]] <- mcmc(lik = con.lik,
                  x.init = runif(min=0, max=1,
                                 n=length(argnames(con.lik))),
                  nsteps = iter,
                  w = w,
                  prior = prior,
                  upper = upper)
  # convert rate parameter for millions of years
  results[[j]] <- results[[j]] / tree.depth
}
stopCluster(cl)
# save
save.image("../results/w.poly.all.matches.MCMC.RData")
