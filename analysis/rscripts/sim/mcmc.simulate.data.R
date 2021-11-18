# load libraries
library(diversitree)
library(chromePlus)
library(doSNOW)

# load data
load("simulate.data.RData")

# define number of clusters for parallel computing
NumberOfClusters <- 2
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)

# MCMC ---------------------
# define parameters for MCMC
iter.temp <- 20
iter <- 100
prior <- make.prior.exponential(r = .5)
results <- vector(mode = "list", length = 100)
# make place holders for results
results.condition.1 <- vector(mode = "list", length = 5)
results.condition.2 <- vector(mode = "list", length = 5)
results.condition.3 <- vector(mode = "list", length = 5)
results.condition.4 <- vector(mode = "list", length = 5)
results.condition.5 <- vector(mode = "list", length = 5)

names(results.condition.1) <- names(trees)
names(results.condition.2) <- names(trees)
names(results.condition.3) <- names(trees)
names(results.condition.4) <- names(trees)
names(results.condition.5) <- c(2,4,6,8,10)

# condition 1 ----
for(i in 1:5){
  x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
    # make the initial data frame
    MCMC.dat <- NULL
    MCMC.dat <-  data.frame(species = names(chrom.traits.1[[i]][[j]]),
                            Chroms = chrom.traits.1[[i]][[j]],
                            bin = bin.traits.1[[i]][[j]],
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
    lik <- make.mkn(tree = trees[[i]][[j]],
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
                            constrain = list(drop.poly=F,
                                             drop.demi=T))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- mcmc(lik = con.lik,
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
    results[[j]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior,
                         lower = rep(0,length(argnames(con.lik))))
    
  }
  results.condition.1[[i]] <- x
}

# condition 1 ----




# condition 2 ----
for(i in 1:5){
  x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
    # make the initial data frame
    MCMC.dat <- NULL
    MCMC.dat <-  data.frame(species = names(chrom.traits.2[[i]][[j]]),
                            Chroms = chrom.traits.2[[i]][[j]],
                            bin = bin.traits.2[[i]][[j]],
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
    lik <- make.mkn(tree = trees[[i]][[j]],
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
                            constrain = list(drop.poly=F,
                                             drop.demi=T))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- mcmc(lik = con.lik,
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
    results[[j]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior,
                         lower = rep(0,length(argnames(con.lik))))
    
  }
  results.condition.2[[i]] <- x
}

# condition 2 ----




# condition 3 ----
for(i in 1:5){
  x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
    # make the initial data frame
    MCMC.dat <- NULL
    MCMC.dat <-  data.frame(species = names(chrom.traits.3[[i]][[j]]),
                            Chroms = chrom.traits.3[[i]][[j]],
                            bin = bin.traits.3[[i]][[j]],
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
    lik <- make.mkn(tree = trees[[i]][[j]],
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
                            constrain = list(drop.poly=F,
                                             drop.demi=T))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- mcmc(lik = con.lik,
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
    results[[j]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior,
                         lower = rep(0,length(argnames(con.lik))))
    
  }
  results.condition.3[[i]] <- x
}

# condition 3 ----




# condition 4 ----
for(i in 1:5){
  x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
    # make the initial data frame
    chroms <- NULL
    MCMC.dat <-  data.frame(species = names(chrom.traits.4[[i]][[j]]),
                            Chroms = chrom.traits.4[[i]][[j]],
                            bin = bin.traits.4[[i]][[j]],
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
    lik <- make.mkn(tree = trees[[i]][[j]],
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
                            constrain = list(drop.poly=F,
                                             drop.demi=T))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- mcmc(lik = con.lik,
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
    results[[j]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior,
                         lower = rep(0,length(argnames(con.lik))))
    
  }
  results.condition.4[[i]] <- x
}

# condition 4 ----




# condition 5 ----
for(i in 1:5){
  x <- foreach(j = 1:100, .verbose = T, .packages = c("ape","diversitree", "chromePlus")) %dopar% {
    # make the initial data frame
    chroms <- NULL
    MCMC.dat <-  data.frame(species = names(chrom.traits.5[[i]][[j]]),
                            Chroms = chrom.traits.5[[i]][[j]],
                            bin = bin.traits.5[[i]][[j]],
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
    lik <- make.mkn(tree = trees[[5]][[j]],
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
                            constrain = list(drop.poly=F,
                                             drop.demi=T))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- mcmc(lik = con.lik,
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
    results[[j]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior,
                         lower = rep(0,length(argnames(con.lik))))
    
  }
  results.condition.5[[i]] <- x
}

# condition 5 ----


# stop the cluster
stopCluster(cl)
# save
save.image("mcmc.simulate.data.RData")

# MCMC ---------------------