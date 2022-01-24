# Terrence Sylvester
# pradakshanas@gmail.com
# 11 January 2022

# load libraries
library(chromePlus)
library(diversitree)
library(phytools)
library(doSNOW)

# define number of clusters for parallel computing
NumberOfClusters <- 2
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)

# set file pars
rr <- 2
# load data name
load.name <- paste("../data/simulate.data.rr-",
                   rr,
                   ".RData",
                   sep = "")
# read in data
load(load.name)
# define parameters for MCMC
iter.temp <- 20
iter <- 100
prior <- make.prior.exponential(r = .5)
ntip <- c(50,100,150,200,250)
# run loop
for(i in 1:6){
  condition <- i
  for(ii in 1:5){
    # read rates
    rates.name <- NULL
    dat <- NULL
    rates.name <- paste("../results/post.burnin/post.burnin.condition.",
                        condition,
                        ".tree",
                        ntip[ii],
                        ".rr.",
                        rr,
                        ".csv",
                        sep = "")
    dat <- read.csv(rates.name, as.is = T)
    for(iii in 1:100){
      # get tree
      tree <- NULL
      tree <- trees[[ii]][[iii]]
      emp.bin.state <- NULL
      if(condition == 1){
        emp.bin.state <- bin.traits.1[[ii]][[iii]]
      }
      if(condition == 2){
        emp.bin.state <- bin.traits.2[[ii]][[iii]]
      }
      if(condition == 3){
        emp.bin.state <- bin.traits.3[[ii]][[iii]]
      }
      if(condition == 4){
        emp.bin.state <- bin.traits.4[[ii]][[iii]]
      }
      if(condition == 5){
        emp.bin.state <- bin.traits.5[[ii]][[iii]]
      }
      if(condition == 6){
        emp.bin.state <- bin.traits.6[[ii]][[iii]]
      }
      ### simulate binary trait data ###
      # get the initial transition matrix
      Q <- matrix(data = 0, 
                  nrow = 2, 
                  ncol = 2)
      # colnames and rownames Q
      colnames(Q) <- rownames(Q) <- c(0,1)
      # Fill in Q
      Q[1,2] <- mean(dat$tran12[((i*50) - 49):(1*50)])
      Q[2,1] <- mean(dat$tran21[((i*50) - 49):(1*50)])
      diag(Q) <- -rowSums(Q)
      # get the ancestral state of the binary character
      root <- fitMk(tree = tree,
                    x = emp.bin.state,
                    fixedQ = Q,
                    pi = "fitzjohn")
      for(iiii in 1:100){ ### multi core step ###
        # simulate binary traits
        sim.bin <- NULL
        sim.bin <- sim.character(tree = tree,
                                 pars = Q,
                                 model = "mkn",
                                 x0 = as.numeric(sample(names(root$pi),
                                                        1, 
                                                        prob = root$pi))+1)
        # change simulated binary characters back to zeros and ones
        sim.bin[sim.bin == 1] <- 0
        sim.bin[sim.bin == 2] <- 1
        
        # get chroms data
        chroms <- NULL
        if(condition == 1){
          chroms <- chrom.traits.1[[ii]][[iii]]
        }
        if(condition == 2){
          chroms <- chrom.traits.2[[ii]][[iii]]
        }
        if(condition == 3){
          chroms <- chrom.traits.3[[ii]][[iii]]
        }
        if(condition == 4){
          chroms <- chrom.traits.4[[ii]][[iii]]
        }
        if(condition == 5){
          chroms <- chrom.traits.5[[ii]][[iii]]
        }
        if(condition == 6){
          chroms <- chrom.traits.6[[ii]][[iii]]
        }
        # make the initial data frame
        MCMC.dat <- NULL
        MCMC.dat <-  data.frame(species = tree$tip.label,
                                chroms = chroms,
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
        
        # store results
        dir.name <- NULL
        dir.name <- paste("../results/nulDist/Ntip",
                          ntip[ii],
                          "/Condition",
                          condition,
                          sep = "")
        if(dir.exists(dir.name) == F){
          dir.create(dir.name, recursive = T)
        }
        # make a name for each file
        save.file.name <- paste("Null.dist.condition.",
                                condition,
                                ".ntip.",
                                ntip[ii],
                                ".tree.",
                                iii,
                                ".sim.",
                                iiii,
                                ".rr-",
                                rr,
                                ".csv",
                                sep = "")
        mcmc(lik = con.lik,
             x.init = runif(min=0, max=1,
                            n=length(argnames(con.lik))),
             nsteps = iter,
             w = w,
             prior = prior,
             lower = rep(0,length(argnames(con.lik))),
             save.file = paste(dir.name,"/",save.file.name, sep = ""),
             save.every = 1)
        x <- NULL
      }
    }
  }
}

