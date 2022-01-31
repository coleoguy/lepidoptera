# Terrence Sylvester
# pradakshanas@gmail.com
# january 25th 2022

# load libraries
library(chromePlus)
library(diversitree)
library(phytools)
library(coda)

#### simulate data ####
# simulate a tree
tree <- tree.bd(pars = c(3, 1), max.taxa = 50)
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# simulate chromosomes and binary trait data
simDat <-  simChrom(tree, 
                    pars=c(.75,     # gains at state 0
                           1*.75,   # gains at state 1
                           .75,     # loss at state 0
                           1*.75,   # loss at state 1
                           0,       # demiploidy at state 0
                           0,       # demiploidy at state 1
                           0,       # polyploidy at state 0
                           0,       # polyploidy at state 1
                           .5,      # transition from state 0 to state 1
                           .5,      # transition from state 1 to state 0
                           10,      # root chromosome number
                           0),      # root state
                    limits = c(1, 100), model = "ChromPlus")
#### simulate data ####


#### run MCMC ####
# define parameters for MCMC
iter.temp <- 20
iter <- 100
prior <- make.prior.exponential(r = .5)

# make data frame
MCMC.dat <- NULL
MCMC.dat <-  data.frame(species = names(simDat$chrom.num),
                        Chroms = simDat$chrom.num,
                        bin = simDat$binary.state,
                        stringsAsFactors = F,
                        row.names = NULL)
# get the range of chromosome number
rng <- NULL
rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
         range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
# convert the data frame to diversitree usable format
chrom.mat <- NULL
chrom.mat <- datatoMatrix(x = MCMC.dat,
                          range = rng,
                          hyper = T)

# make the likelihood function
lik <- NULL
lik <- make.mkn(tree = tree,
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = F,
                control = list(method="ode"))
# constrain the likelihood function
con.lik <- NULL
con.lik <- constrainMkn(data = chrom.mat,
                        lik = lik,
                        hyper = T,
                        polyploidy = F,
                        verbose = F,
                        oneway = F,
                        constrain = list(drop.poly=T,
                                         drop.demi=T))

# run the initial MCMC to get parameter values for w
temp <- NULL
temp <- diversitree::mcmc(lik = con.lik,
                          x.init = runif(min=0, max=1,
                                         n=length(argnames(con.lik))),
                          prior = prior,
                          nsteps = iter.temp,
                          w = 1,
                          lower = rep(0,length(argnames(con.lik))))
# get values for w
w <- NULL
w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
# run MCMC
results <- NULL
results <- diversitree::mcmc(lik = con.lik,
                             x.init = runif(min=0, max=1,
                                            n=length(argnames(con.lik))),
                             nsteps = iter,
                             w = w,
                             prior = prior,
                             lower = rep(0,length(argnames(con.lik))))
#### run MCMC ####

#### process MCMC ####
# get post burnin
plot(results$p, type = "l")
postBurnResults <- results[c(51:100),]
#### process MCMC ####


#### emperical p value ####
resultsP <- vector(mode = "list", length = 100)
# get the transistion matrix
tmat <- matrix(data = 0, nrow = 2, ncol = 2)
colnames(tmat) <- rownames(tmat) <- c(0,1)
# fill tmat
tmat[1,2] <- mean(postBurnResults$tran12)
tmat[2,1] <- mean(postBurnResults$tran21)
diag(tmat) <- -rowSums(tmat)
# get the root probabilities of the binary trait
root <- fitMk(tree = tree, 
              x = simDat$binary.state,
              fixedQ = tmat,
              pi = "fitzjohn")
for(i in 1:100){
  # simulate binary traits
  sim.bin <- NULL
  sim.bin <- sim.character(tree = tree,
                           pars = tmat,
                           model = "mkn",
                           x0 = as.numeric(sample(names(root$pi),
                                                  1, 
                                                  prob = root$pi))+1)
  # change simulated binary characters back to zeros and ones
  sim.bin[sim.bin == 1] <- 0
  sim.bin[sim.bin == 2] <- 1
  
  # make the initial data frame
  MCMC.dat <- NULL
  MCMC.dat <-  data.frame(species = tree$tip.label,
                          chroms = simDat$chrom.num,
                          bin = sim.bin,
                          stringsAsFactors = F,
                          row.names = NULL)
  # get the range of chromosome number
  rng <- NULL
  rng <- c(range(MCMC.dat$chroms, na.rm = T)[1] - 1,
           range(MCMC.dat$chroms, na.rm = T)[2] + 1)
  # convert the data frame to diversitree usable format
  chrom.mat <- NULL
  chrom.mat <- datatoMatrix(x = MCMC.dat,
                            range = rng,
                            hyper = T)
  
  # make the likelihood function
  lik <- NULL
  lik <- make.mkn(tree = tree,
                  states = chrom.mat,
                  k = ncol(chrom.mat),
                  strict = F,
                  control = list(method="ode"))
  # constrain the likelihood function
  con.lik <- NULL
  con.lik <- constrainMkn(data = chrom.mat,
                          lik = lik,
                          hyper = T,
                          polyploidy = F,
                          verbose = F,
                          oneway = F,
                          constrain = list(drop.poly=T,
                                           drop.demi=T))
  
  # run the initial MCMC to get parameter values for w
  temp <- NULL
  temp <- diversitree::mcmc(lik = con.lik,
                            x.init = runif(min=0, max=1,
                                           n=length(argnames(con.lik))),
                            prior = prior,
                            nsteps = iter.temp,
                            w = 1,
                            lower = rep(0,length(argnames(con.lik))))
  # get values for w
  w <- NULL
  w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
  # run MCMC
  resultsP[[i]] <-  diversitree::mcmc(lik = con.lik,
                                      x.init = runif(min=0, max=1,
                                                     n=length(argnames(con.lik))),
                                      nsteps = iter,
                                      w = w,
                                      prior = prior,
                                      lower = rep(0,length(argnames(con.lik))))
  
}

#### process MCMC ####
# get postburnin
y <- c()
y <- resultsP[[1]][51:100,]
for(i in 2:100){
  y <- rbind(y,resultsP[[i]][51:100,])
}
postBurnResultsP <- y
#### process MCMC ####

#### plotting ####
# make a table to get the HPD intervals for each MCMC we did
HPD.tab.fusions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
HPD.tab.fissions <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
# HPD.tab.poly <- as.data.frame(matrix(data = NA, nrow = 10, ncol = 4))
# get column names
colnames(HPD.tab.fusions) <- colnames(HPD.tab.fissions) <- c("tree", "HPD-low", "HPD-high", "pass-zero")
# make a table to get the distribution of mean delta r statistic
Delta.R <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 4))
colnames(Delta.R) <- c("Sim", "asc", "desc", "pol")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# fill in data
for(i in 1:100){
  # get the working tree
  HPD.tab.fusions$tree[i] <- i
  HPD.tab.fissions$tree[i] <- i
  # HPD.tab.poly$tree[i] <- i
  # get the lower limit of HPD
  HPD.tab.fusions$`HPD-low`[i] <- HPDinterval(as.mcmc(postBurnResultsP$asc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$asc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  HPD.tab.fissions$`HPD-low`[i] <- HPDinterval(as.mcmc(postBurnResultsP$desc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$desc2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  # HPD.tab.poly$`HPD-low`[i] <- HPDinterval(as.mcmc(postBurnResultsP$pol1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$pol2[pb.seq[i]:(pb.seq[i]+49)]))[1]
  # get the higher limit of HPD
  HPD.tab.fusions$`HPD-high`[i] <- HPDinterval(as.mcmc(postBurnResultsP$asc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$asc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  HPD.tab.fissions$`HPD-high`[i] <- HPDinterval(as.mcmc(postBurnResultsP$desc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$desc2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  # HPD.tab.poly$`HPD-high`[i] <- HPDinterval(as.mcmc(postBurnResultsP$pol1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$pol2[pb.seq[i]:(pb.seq[i]+49)]))[2]
  # get the mean of the delta R statistic for each run
  Delta.R$Sim[i] <- i
  Delta.R$asc[i] <- mean(abs(postBurnResultsP$asc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$asc2[pb.seq[i]:(pb.seq[i]+49)]))
  Delta.R$desc[i] <- mean(abs(postBurnResultsP$desc1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$desc2[pb.seq[i]:(pb.seq[i]+49)]))
  # Delta.R$pol[i] <- mean(abs(postBurnResultsP$pol1[pb.seq[i]:(pb.seq[i]+49)] - postBurnResultsP$pol2[pb.seq[i]:(pb.seq[i]+49)]))
  
  if(HPD.tab.fusions$`HPD-low`[i] < 0 & HPD.tab.fusions$`HPD-high`[i] > 0){
    HPD.tab.fusions$`pass-zero`[i] <- T  
  }else{
    HPD.tab.fusions$`pass-zero`[i] <- F
  }
  if(HPD.tab.fissions$`HPD-low`[i] < 0 & HPD.tab.fissions$`HPD-high`[i] > 0){
    HPD.tab.fissions$`pass-zero`[i] <- T  
  }else{
    HPD.tab.fissions$`pass-zero`[i] <- F
  }
  # if(HPD.tab.poly$`HPD-low`[i] < 0 & HPD.tab.poly$`HPD-high`[i] > 0){
  #   HPD.tab.poly$`pass-zero`[i] <- T  
  # }else{
  #   HPD.tab.poly$`pass-zero`[i] <- F
  # }
}
sum(HPD.tab.fusions$`pass-zero`)
sum(HPD.tab.fissions$`pass-zero`)
# sum(HPD.tab.poly$`pass-zero`)



# plot the postBurnResults
par(mfcol = c(1,2))
# fission
hist(Delta.R$asc, 
     xlab = expression(paste("| ",Delta, "R"[fission]," |")),
     main = "",
     breaks = 25,
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(postBurnResults$asc1 - postBurnResults$asc2)),
       col = "red",
       lwd = 2)
mtext(text = "A",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# text(x = min(density(Delta.R$asc)$x),
#      y = max(density(Delta.R$asc)$y),
#      label = paste("p-value",sum(Delta.R$asc > mean(abs(postBurnResults$asc1 - postBurnResults$asc2))) / 100),
#      pos = 4,
#      cex = 1.2) 
# fusion
hist(Delta.R$desc, 
     xlab = expression(paste("| ",Delta, "R"[fusion]," |")),
     main = "",
     breaks = 25,
     cex.lab = 1.5,
     lwd = 2)
abline(v = mean(abs(postBurnResults$desc1 - postBurnResults$desc2)),
       col = "red",
       lwd = 2)
mtext(text = "B",
      line = 0,
      outer = F,
      side = 3,
      adj = 0,
      cex = 1)
# text(x = max(density(Delta.R$desc)$x),
#      y = max(density(Delta.R$desc)$y),
#      label = paste("p-value",sum(Delta.R$desc > mean(abs(postBurnResults$desc1 - postBurnResults$desc2))) / 100),
#      pos = 2,
#      cex = 1.2) 

#### plotting ####
