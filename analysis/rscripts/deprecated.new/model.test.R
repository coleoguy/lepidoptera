# load libraries
library(ape)
library(phytools)
library(chromePlus)
library(diversitree)
library(doSNOW)

# define number of clusters for parrallel computing
NumberOfCluster <- 4
cl <- makeCluster(NumberOfCluster, outfile="")
registerDoSNOW(cl)

# parameters for MCMC
iter <- 100
iter.temp <- 20
prior <- make.prior.exponential(r=2)
upper <- 50
results <- vector(mode = "list", length = 10)

# get helper functions
source("helper.functions.R")

# read in post burnin portion of the mcmc 
post.burnin <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
# mean genaralist to specialist rate
q01 <- mean(post.burnin$tran21) # in mcmc state 1 was specialist and 2 was generalist
# mean specialist to generalist rate
q10 <- mean(post.burnin$tran12)

# run the MCMC
x <- foreach(j = 1:100) %dopar% {
  library(ape)
  library(phytools)
  library(chromePlus)
  library(diversitree)
  # read in data
  dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
  trees <- read.tree("../data/trees/processed.trees.new")
  # combine genus and species names to get the full name in chroms dataset
  # fill in species for data
  dat$SpecisName <- paste(dat$other.names.genera,
                          dat$other.names.species,
                          sep = "_")
  dat$SpecisName[dat$SpecisName == "_"] <- paste(dat$Genus[dat$SpecisName == "_"],
                                                 dat$species[dat$SpecisName == "_"],
                                                 sep = "_")
  # randomly sample a single tree for tha analysis
  sampled.tree <- trees[[j]]
  # sample a single chromosome number when there is a range
  new.dat <- chromSampler(dat = dat)
  # remove species that lack male chromosome number data
  new.dat <-  new.dat[!(is.na(new.dat$male2N)),]
  # get the overlap between the traits and the phylogeny
  overlap <- sp.matches(dat = new.dat, trees = trees,use.sub.species = F)
  # proccessed dat
  processed.dat <- overlap$chroms
  # rename tree and hosts species names
  tree <- keep.tip(sampled.tree, overlap$name_corrections$name_on_tree)
  # rename tree tip and host species name to reflect genera level matches
  for(i in 1:Ntip(tree)){
    tree$tip.label[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == tree$tip.label[i]]
  }
  # get the tip states of a neutrally evolving trait
  tip.states <- sim.character(tree = tree,
                              pars = c(q01,q10), #(0 to 1, 1 to 0)
                              model = "mk2", 
                              x0 = sample(c(0,1),1)) # root state
  # make a new data table which has species names chromosome number and host type
  MCMC.dat <- as.data.frame(matrix(data = NA,
                                   nrow = nrow(processed.dat),
                                   ncol = 3))
  colnames(MCMC.dat) <- c("Species", "Chroms", "State")
  MCMC.dat$Species <- processed.dat$SpecisName
  MCMC.dat$Chroms <- processed.dat$male2N / 2
  for(i in 1:nrow(MCMC.dat)){
    MCMC.dat$State[i] <- tip.states[names(tip.states) == MCMC.dat$Species[i]]
  }
  # tree is not utrametric (eventhough it is a BEAST output). Make that correction
  tree <- force.ultrametric(tree)
  # remove unwanted tips
  tree <- drop.tip(tree, processed.dat$SpecisName[is.na(processed.dat$male2N)])
  # make tree unit length
  tree.depth <- max(branching.times(tree))
  tree$edge.length <- tree$edge.length / tree.depth
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
  # get values for w
  temp <- mcmc(lik = con.lik,
               x.init = runif(min=0, max=1,
                              n=length(argnames(con.lik))),
               nsteps = iter.temp,
               w = 1,
               prior = prior,
               upper = upper)
  w <- diff(sapply(temp[11:20, 2:9], quantile, c(.05, .95)))
  # run the mcmc
  results[[j]] <- mcmc(lik = con.lik,
                       x.init = runif(min=0, max=1,
                                      n=length(argnames(con.lik))),
                       nsteps = iter,
                       w = w,
                       prior = prior,
                       upper = upper)
  results[[j]] <- results[[j]] / tree.depth
}

stopCluster(cl)
# save
save.image("../results/model.test.RData")
