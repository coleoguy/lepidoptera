# Terrence Sylvester
# 29 June 2021
# pradakshanas@gmail.com
# model adequacy testing

# load libraries
library(ape)
library(phytools)
library(chromePlus)
library(diversitree)

# get helper functions
source("helper.functions.R")

# read in data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
# combine genus and species names to get the full name in hosts dataset
host$binomial <- paste(host$genus, host$species, sep = "_")

# read in post burnin portion of the mcmc 
post.burnin <- read.csv("../results/w.poly.all.matches.post.burnin.csv", as.is = T)
# mean genaralist to specialist rate
q01 <- mean(post.burnin$tran21) # in mcmc state 1 was specialist and 2 was generalist
# mean specialist to generalist rate
q10 <- mean(post.burnin$tran12)
# mean rate of chromosome fission at state 1 
state.01.asc <- mean(post.burnin$asc1)
# mean rate of chromosom fission at state 2
state.02.asc <- mean(post.burnin$asc2)
# mean rate of chromosome fusion at state 1
state.01.desc <- mean(post.burnin$desc1)
# mean rate of chromosome fusuin at state 2
state.02.desc <- mean(post.burnin$desc2)
# mean rate of polyploidy at state 1
state.01.pol <- mean(post.burnin$pol1)
# mean rate of polyploidy at state 2
state.02.pol <- mean(post.burnin$pol2)

# place holder for final results
results.final <- vector(mode = "list", length = 10)
names(results.final) <- paste("tree", 1:10, sep = "")
# do the simulations
for(j in 1:10){
  # place holder for results
  results.chroms <- vector(mode = "list", length = 1000)
  results.states <- vector(mode = "list", length = 1000)
  pic.chroms <- vector(mode = "list", length = 1000)
  pic.states <- vector(mode = "list", length = 1000)
  results.temp <- c()
  names(results.chroms) <- names(results.states) <- names(pic.chroms) <- names(pic.states) <- paste("sim", 1:1000, sep = "")
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
  host.sub <- host[host$binomial %in% overlap$name_corrections$name_on_tree,]
  for(i in 1:nrow(host.sub)){
    host.sub$binomial[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == host.sub$binomial[i]]
  }
  # get named states
  f.type <- host.sub$type
  names(f.type) <- host.sub$binomial
  
  # make a new data table which has species names chromosome number and host type
  MCMC.dat <- as.data.frame(matrix(data = NA,
                                   nrow = nrow(processed.dat),
                                   ncol = 4))
  colnames(MCMC.dat) <- c("Species", "Chroms", "HostType","State")
  MCMC.dat$Species <- processed.dat$binomial
  MCMC.dat$Chroms <- processed.dat$male2N / 2
  for(i in 1:nrow(MCMC.dat)){
    MCMC.dat$HostType[i] <- host.sub$type[host.sub$binomial == MCMC.dat$Species[i]]
  }
  
  # tree is not utrametric (eventhough it is a BEAST output). Make that correction
  tree <- force.ultrametric(tree)
  # remove unwanted tips
  tree <- drop.tip(tree, processed.dat$SpecisName[is.na(processed.dat$male2N)])
  # get tree depth
  tree.depth <- max(branching.times(tree))
  # scale tree to unit length
  unit.tree <-tree 
  unit.tree$edge.length / tree.depth
  # properly name each state in the qmatrix
  for(i in 1:nrow(MCMC.dat)){
    if (max(MCMC.dat$Chroms) < 100) 
      pad <- 2
    if (max(MCMC.dat$Chroms) >= 100) 
      pad <- 3
    if (max(MCMC.dat$Chroms) < 10) 
      pad <- 1
    if(MCMC.dat$HostType[i] == 1){
      MCMC.dat$State[i] <-  sprintf(paste("%0", pad, "d", 
                                          sep = ""), which(sort(min(MCMC.dat$Chroms):max(MCMC.dat$Chroms)) == MCMC.dat$Chroms[i]) + length(min(MCMC.dat$Chroms):max(MCMC.dat$Chroms)))
    }else{
      MCMC.dat$State[i] <-  sprintf(paste("%0", pad, "d", 
                                          sep = ""), which(sort(min(MCMC.dat$Chroms):max(MCMC.dat$Chroms)) == MCMC.dat$Chroms[i]))
    }
  }
  # buld qmatrix
  # get pars
  pars = c(state.01.asc,   # ascending at state 1
           state.02.asc,   # ascending at state 2
           state.01.desc,  # descending at state 1
           state.02.desc,  # descending at state 2
           0,              # demiploidy at state 1
           0,              # demiploidy at state 2
           state.01.pol,   # polyploidy at state 1
           state.02.pol,   # polyploidy at state 2
           q10,            # transision from state 1 to state 2
           q01)            # transision from state 2 to state 1
  q <- buildQ(pars = pars, limits = range(MCMC.dat$Chroms))
  # get the root of the tree
  for(i in 1:nrow(q)){
    q[i,i] <- -sum(q[i,])
  }
  # make two copies of the qmatrix to be used in two different functions
  q.model.fit <- q
  q.sim.char <- q  
  # prepeare inputs for fitMK
  chromMat <-  datatoMatrix(x = MCMC.dat, range = range(MCMC.dat$Chroms),hyper = T)
  colnames(q.model.fit) <- rownames(q.model.fit) <- colnames(chromMat)
  # get root states
  model.fit <- fitMk(tree = tree, x = chromMat,fixedQ = q.model.fit, pi="fitzjohn")
  # rename the colnames of q matrix used in sumulating datasets
  colnames(q.sim.char) <- rownames(q.sim.char) <- 1:ncol(q.sim.char)
  names(model.fit$pi) <- 1:ncol(q.sim.char)
  # get the pics of observed data
  named.chroms <- MCMC.dat$Chroms
  named.states <- MCMC.dat$HostType
  names(named.chroms) <- names(named.states) <- MCMC.dat$Species
  pic.chroms.obs <- pic(x = named.chroms, phy = unit.tree,var.contrasts = T)
  pic.states.obs <- pic(x = named.states, phy = unit.tree, var.contrasts = T)
  # simulate datasets
  for(k in 1:1000){
    print(paste("Tree", j, "simulation", k))
    # simulate datasets
    sim.chroms <- sim.character(tree = tree,
                                pars = q.sim.char,
                                x0 = as.numeric(sample(names(model.fit$pi),size = 1,prob = model.fit$pi)),
                                model = "mkn")
    # get states
    states <- vector(mode = "numeric", length = length(sim.chroms))
    states[sim.chroms > max(MCMC.dat$Chroms)] <- 0
    states[sim.chroms <= max(MCMC.dat$Chroms)] <- 1
    # get chroms
    sim.chrom.states <- vector(mode = "numeric", length = length(sim.chroms))
    sim.chrom.states[sim.chroms > max(MCMC.dat$Chroms)] <- sim.chroms[sim.chroms > max(MCMC.dat$Chroms)] - length(min(MCMC.dat$Chroms):max(MCMC.dat$Chroms)) + (min(MCMC.dat$Chroms) - 1)
    sim.chrom.states[sim.chroms <= max(MCMC.dat$Chroms)] <- sim.chroms[sim.chroms <= max(MCMC.dat$Chroms)] + (min(MCMC.dat$Chroms) - 1)
    results.chroms[[k]] <- sim.chrom.states
    results.states[[k]] <- states
    # pic
    pic.states[[k]] <- pic(x = states, phy = unit.tree,var.contrasts = T)
    pic.chroms[[k]] <- pic(x = sim.chroms, phy = unit.tree,var.contrasts = T)
  }
  results.temp <- list(pic.chroms.obs,
                       pic.states.obs,
                       results.chroms,
                       results.states,
                       pic.chroms,
                       pic.states)
  names(results.temp) <- c("pic.chroms.obs",
                           "pic.states.obs",
                           "chroms",
                           "states",
                           "pic.sim.chroms",
                           "pic.sim.states")
  results.final[[j]] <- results.temp
}


save.image("../results/model.adequacy.testing.RData")

# analysis and plotting
# get the statistics
# stats of the obs data
# chroms
Msig.obs.chrom <- vector(mode = "numeric", length = 10)
Cvar.obs.chrom <- vector(mode = "numeric", length = 10)
Svar.obs.chrom <- vector(mode = "numeric", length = 10)
Dcdf.obs.chrom <- vector(mode = "numeric", length = 10)
# state
Msig.obs.state <- vector(mode = "numeric", length = 10)
Cvar.obs.state <- vector(mode = "numeric", length = 10)
Svar.obs.state <- vector(mode = "numeric", length = 10)
Dcdf.obs.state <- vector(mode = "numeric", length = 10)
# stats of the sim data
# chroms
Msig.sim.chrom <- vector(mode = "numeric", length = 10000)
Cvar.sim.chrom <- vector(mode = "numeric", length = 10000)
Svar.sim.chrom <- vector(mode = "numeric", length = 10000)
Dcdf.sim.chrom <- vector(mode = "numeric", length = 10000)
# state
Msig.sim.state <- vector(mode = "numeric", length = 10000)
Cvar.sim.state <- vector(mode = "numeric", length = 10000)
Svar.sim.state <- vector(mode = "numeric", length = 10000)
Dcdf.sim.state <- vector(mode = "numeric", length = 10000)
# get stats
counter <- 1
for(i in 1:10){
  # chroms
  Msig.obs.chrom[i] <- mean(results.final[[i]][[1]][,1] ^ 2)
  Cvar.obs.chrom[i] <- sd(abs(results.final[[i]][[1]][,1])) / mean(abs(results.final[[i]][[1]][,1]))
  Svar.obs.chrom[i] <- lm(abs(results.final[[i]][[1]][,1]) ~ results.final[[i]][[1]][,2])[[1]][2]
  Dcdf.obs.chrom[i] <- ks.test(x = results.final[[i]][[1]][,1],
                               y = rnorm(n = length(results.final[[i]][[1]][,1]),
                                         sd = sqrt(mean(results.final[[i]][[1]][,1]^ 2))))$statistic
  # state
  Msig.obs.state[i] <- mean(results.final[[i]][[2]][,1] ^ 2)
  Cvar.obs.state[i] <- sd(abs(results.final[[i]][[2]][,1])) / mean(abs(results.final[[i]][[2]][,1]))
  Svar.obs.state[i] <- lm(abs(results.final[[i]][[2]][,1]) ~ results.final[[i]][[2]][,2])[[1]][2]
  Dcdf.obs.state[i] <- ks.test(x = results.final[[i]][[2]][,1],
                               y = rnorm(n = length(results.final[[i]][[2]][,1]),
                                         sd = sqrt(mean(results.final[[i]][[2]][,1]^ 2))))$statistic
  for(j in 1:1000){
    # chroms
    Msig.sim.chrom[counter] <- mean(results.final[[i]][[5]][[j]][,1] ^ 2)
    Cvar.sim.chrom[counter] <- sd(abs(results.final[[i]][[5]][[j]][,1])) / mean(abs(results.final[[i]][[5]][[j]][,1]))
    Svar.sim.chrom[counter] <- lm(abs(results.final[[i]][[5]][[j]][,1]) ~ results.final[[i]][[5]][[j]][,2])[[1]][2]
    Dcdf.sim.chrom[counter] <- ks.test(x = results.final[[i]][[5]][[j]][,1],
                                       y = rnorm(n = length(results.final[[i]][[5]][[j]][,1]),
                                                 sd = sqrt(mean(results.final[[i]][[5]][[j]][,1]^ 2))))$statistic
    # state
    Msig.sim.state[counter] <- mean(results.final[[i]][[6]][[j]][,1] ^ 2)
    Cvar.sim.state[counter] <- sd(abs(results.final[[i]][[6]][[j]][,1])) / mean(abs(results.final[[i]][[6]][[j]][,1]))
    Svar.sim.state[counter] <- lm(abs(results.final[[i]][[6]][[j]][,1]) ~ results.final[[i]][[6]][[j]][,2])[[1]][2]
    Dcdf.sim.state[counter] <- ks.test(x = results.final[[i]][[6]][[j]][,1],
                                       y = rnorm(n = length(results.final[[i]][[6]][[j]][,1]),
                                                 sd = sqrt(mean(results.final[[i]][[6]][[j]][,1]^ 2))))$statistic
    counter <- counter + 1
  }
}

# plot stats
par(mfrow = c(2,4))
# Msig chrom
plot(density(Msig.sim.chrom),
     main = "",
     xlab = "Msig", 
     xlim = c(min(c(Msig.obs.chrom, density(Msig.sim.chrom)$x)),
              max(c(Msig.obs.chrom, density(Msig.sim.chrom)$x))),
     lwd = 2)
abline(v = sample(Msig.obs.chrom,1), col = "red", lty = 2,
       lwd = 2)
# Cvar chrom
plot(density(Cvar.sim.chrom),
     main = "",
     xlab = "Cvar", 
     xlim = c(min(c(Cvar.obs.chrom, density(Cvar.sim.chrom)$x)),
                             max(c(Cvar.obs.chrom, density(Cvar.sim.chrom)$x))),
     lwd = 2)
abline(v = sample(Cvar.obs.chrom,1), col = "red", lty = 2,
       lwd = 2)
# Svar chrom
plot(density(Svar.sim.chrom),
     main = "",
     xlab = "Svar", 
     xlim = c(min(c(Svar.obs.chrom, density(Svar.sim.chrom)$x)),
              max(c(Svar.obs.chrom, density(Svar.sim.chrom)$x))),
     lwd = 2)
abline(v = sample(Svar.obs.chrom,1), col = "red", lty = 2,
       lwd = 2)
# Dcdf chrom
plot(density(Dcdf.sim.chrom),
     main = "",
     xlab = "Msig", 
     xlim = c(min(c(Dcdf.obs.chrom, density(Dcdf.sim.chrom)$x)),
              max(c(Dcdf.obs.chrom, density(Dcdf.sim.chrom)$x))),
     lwd = 2)
abline(v = sample(Dcdf.obs.chrom,1), col = "red", lty = 2,
       lwd = 2)

# states
# Msig state
plot(density(Msig.sim.state),
     main = "",
     xlab = "Msig", 
     xlim = c(min(c(Msig.obs.state, density(Msig.sim.state)$x)),
              max(c(Msig.obs.state, density(Msig.sim.state)$x))),
     lwd = 2)
abline(v = sample(Msig.obs.state,1), col = "red", lty = 2,
       lwd = 2)
# Cvar state
plot(density(Cvar.sim.state),
     main = "",
     xlab = "Cvar", 
     xlim = c(min(c(Cvar.obs.state, density(Cvar.sim.state)$x)),
              max(c(Cvar.obs.state, density(Cvar.sim.state)$x))),
     lwd = 2)
abline(v = sample(Cvar.obs.state,1), col = "red", lty = 2,
       lwd = 2)
# Svar state
plot(density(Svar.sim.state),
     main = "",
     xlab = "Svar", 
     xlim = c(min(c(Svar.obs.state, density(Svar.sim.state)$x)),
              max(c(Svar.obs.state, density(Svar.sim.state)$x))),
     lwd = 2)
abline(v = sample(Svar.obs.state,1), col = "red", lty = 2,
       lwd = 2)
# Dcdf state
plot(density(Dcdf.sim.state),
     main = "",
     xlab = "Msig", 
     xlim = c(min(c(Dcdf.obs.state, density(Dcdf.sim.state)$x)),
              max(c(Dcdf.obs.state, density(Dcdf.sim.state)$x))),
     lwd = 2)
abline(v = sample(Dcdf.obs.state,1), col = "red", lty = 2,
       lwd = 2)

par(mfrow = c(2,3))
# proportion of specialists in emperical dataset
obs.prop.div <-  sum(MCMC.dat$HostType)/length(MCMC.dat$HostType)
# proportion of specialists in simulated datasets
sim.prop.div <- c()
for(i in 1:10){
  for(j in 1:1000){
    sim.prop.div <-  c(sim.prop.div, sum(results.final[[i]]$states[[j]])/length(results.final[[i]]$states[[j]]))
  }
}
plot(density(sim.prop.div),
     xlab = "Proportion of specialist species",
     main = "",
     lwd = 2)
abline(v = obs.prop.div,
       col = "red",
       lwd = 2)
# observed coefficient of variance in chromosome number
obs.coef.var <- sd(MCMC.dat$Chroms) / mean(MCMC.dat$Chroms)
# simulated coefficient of variance in chromosome number
sim.coef.var <- c()
for(i in 1:10){
  for(j in 1:1000){
    sim.coef.var <- c(sim.coef.var,(sd(results.final[[i]]$chroms[[j]])/mean(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.coef.var),
     xlab = "Coefficient of variance in chromosome number",
     main = "",      lwd = 2)
abline(v = obs.coef.var,
       col = "red",
       lwd = 2)
# observed variance in chromosome number
obs.var <- var(MCMC.dat$Chroms)
# simulated variance in chromosome number
sim.var <- c()
for(i in 1:10){
  for(j in 1:1000){
    sim.var <- c(sim.var,(var(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.var),
     xlab = "Variance in chromosome number",
     main = "",
     lwd = 2)
abline(v = obs.var,
       col = "red",
       lwd = 2)
# observed SD of chromosome number
obs.sd <- sd(MCMC.dat$Chroms)
# simulated SD of chromosome number
sim.sd <- c()
for(i in 1:10){
  for(j in 1:1000){
    sim.sd <- c(sim.sd,(sd(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.sd),
     xlab = "Standard diviation in chromosome number",
     main = "",
     lwd = 2)
abline(v = obs.sd,
       col = "red",
       lwd = 2)
# observed mean of chromosome number
obs.mean <- mean(MCMC.dat$Chroms)
# simulated mean of chromosome number
sim.mean <- c()
for(i in 1:10){
  for(j in 1:1000){
    sim.mean <- c(sim.mean,(mean(results.final[[i]]$chroms[[j]])))
  }
}
plot(density(sim.mean),
     xlab = "Mean in chromosome number",
     main = "",
     lwd = 2)
abline(v = obs.mean,
       col = "red",
       lwd = 2)
# get the density distribution of chromosome number
# plot(x = NULL, y = NULL,
#      xlim = c(0,100),
#      ylim = c(0,.2),
#      lwd = 2)
# for(i in 1:10){
#   for(j in 1:1000){
#     lines(density(results.final[[i]]$chroms[[j]], bw = 1),
#           col = rgb(1,0,0,.2),
#           lwd = .2)
#   }
# }
# lines(density(MCMC.dat$Chroms, bw = 1),
#       lwd = 2)






