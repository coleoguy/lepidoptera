# this function will simulate a set of birth death trees ----
simTrees <- function(ntaxa = NULL,
                     ntrees = NULL,
                     birth = NULL,
                     death = NULL){
  trees <- list()
  for(i in 1:length(ntaxa)){
    tree.set <- list()
    for(j in 1:ntrees){
      tree <- NULL
      while(is.null(tree)){
        tree <- tree.bd(pars = c(birth, death), max.taxa = ntaxa[i])
      }
      tree$edge.length <- tree$edge.length/max(branching.times(tree))
      tree.set[[j]] <- tree
    }
    names(tree.set) <- paste("tree", 1:ntrees, sep = "")
    trees[[i]] <- tree.set 
  }
  # name each set in trees
  names(trees) <- paste("nTips",ntaxa, sep = "") 
  return(trees)
}
# This function will perform the ChromPlus analysis on empirical data ----
runMCMC <- function(tree = NULL,
                    chroms = NULL,
                    binary = NULL,
                    args.lik = list(control = "ode",
                                    strict = F),
                    args.conlik = list(hyper = T,
                                       polyploidy = F,
                                       verbose = F,
                                       oneway = F,
                                       drop.poly=T,
                                       drop.demi=T,
                                       symmetric=F,
                                       nometa=F,
                                       meta="ARD"),
                    args.MCMC = list(iter.temp = 20,
                                     iter = 100,
                                     prior = make.prior.exponential(r = .5),
                                     print.every=50)){
  # run time
  start_time <- as.numeric(Sys.time())
  # make initial data frame
  MCMC.dat <- NULL
  MCMC.dat <-  data.frame(species = tree$tip.label,
                          Chroms = chroms,
                          bin = binary,
                          stringsAsFactors = F,
                          row.names = NULL)
  # get the range of chromosome number
  rng <- NULL
  rng <- c(range(MCMC.dat$Chroms, na.rm = T)[1] - 1,
           range(MCMC.dat$Chroms, na.rm = T)[2] + 1)
  # convert the data frame to diversitree usable format
  # convert the data frame to diversitree usable format
  chrom.mat <- NULL
  chrom.mat <- chromePlus::datatoMatrix(x = MCMC.dat,
                                        range = rng,
                                        hyper = args.conlik$hyper)
  # make the likelihood function
  lik <- NULL
  lik <- diversitree::make.mkn(tree = tree,
                               states = chrom.mat,
                               k = ncol(chrom.mat),
                               strict = F,
                               control = list(method=args.lik$control))
  # constrain the likelihood function
  con.lik <- NULL
  con.lik <- chromePlus::constrainMkn(data = chrom.mat,
                                      lik = lik,
                                      hyper = args.conlik$hyper,
                                      polyploidy = args.conlik$polyploidy,
                                      verbose = args.conlik$verbose,
                                      oneway = args.conlik$oneway,
                                      constrain = list(drop.poly=args.conlik$drop.poly,
                                                       drop.demi=args.conlik$drop.demi,
                                                       symmetric = args.conlik$symmetric,
                                                       nometa = args.conlik$nometa,
                                                       meta = args.conlik$meta))
  # run the initial MCMC to get parameter values for w
  temp <- NULL
  temp <- diversitree::mcmc(lik = con.lik,
                            x.init = runif(min=0, max=1,
                                           n=length(argnames(con.lik))),
                            prior = args.MCMC$prior,
                            nsteps = args.MCMC$iter.temp,
                            w = 1,
                            lower = rep(0,length(argnames(con.lik))),
                            print.every=5)
  # get values for w
  w <- NULL
  w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
  # run MCMC
  results <- NULL
  results <- diversitree::mcmc(lik = con.lik,
                               x.init = runif(min=0, max=1,
                                              n=length(argnames(con.lik))),
                               nsteps = args.MCMC$iter,
                               w = w,
                               prior = args.MCMC$prior,
                               lower = rep(0,length(argnames(con.lik))),
                               print.every=args.MCMC$print.every)
  end_time <- as.numeric(Sys.time())
  print(paste("run time MCMC =",  round(end_time - start_time,0), "seconds"))
  return(results)
}

# This function will plot the likelihood of the MCMC ----
plotlikMCMC <- function(data, burn = NULL){
  if(is.data.frame(data) == T){
    lik <- which(colnames(data) == "p")
    index <- which(colnames(data) == "i")
    rng <- range(data[,lik])
    plot(y = data[,lik],
         x = data[,index],
         type = "l",
         ylim = rng,
         xlab = "Generation",
         ylab = "Likelihood")
    abline(v = nrow(data) * burn, col = "red")
  }else{
    nElements <- length(data)
    rngMin <- rngMax <- vector(mode = "numeric", length = nElements)
    for(i in 1:nElements){
      rngMin[i] <- range(data[[i]]$p)[1]
      rngMax[i] <- range(data[[i]]$p)[2]
    }
    rng <- range(c(rngMin, rngMax))
    plot(y = data[[1]]$p,
         x = data[[1]]$i,
         type = "l",
         ylim = rng,
         col = viridis::viridis(n = nElements)[1],
         xlab = "Generation",
         ylab = "Likelihood")
    for(i in 2:nElements){
      lines(y = data[[i]]$p,
            x = data[[i]]$i,
            col = viridis::viridis(n = nElements)[i])
    }
    abline(v = nrow(data[[1]]) * burn, col = "red")
  }
}

# This function will get the post burnin from data ----
getPostBurnin <- function(data, burn = NULL){
  if(is.data.frame(data) == T){
    pbrn <- data[((nrow(data)*burn)+1):nrow(data),]
  }else{
    y <- c()
    y <- data[[1]][((nrow(data[[1]])*burn)+1):nrow(data[[1]]),]
    for(i in 2:length(data)){
      y <- rbind(y,data[[i]][((nrow(data[[1]])*burn)+1):nrow(data[[1]]),])
    }
    pbrn <- y  
  }
  return(pbrn)
}

# This function will calculate empirical p value ----
empiricalPcalc <- function(empPostburnin = NULL,
                           simMCMC = NULL,
                           polyploidy = NULL,
                           iter = NULL,
                           nsim = NULL,
                           burn = NULL,
                           plot = F){
  # get mean rate of empirical data
  empMeanASC <- mean(empPostburnin$asc1 - empPostburnin$asc2)
  empMeanDESC <-mean(empPostburnin$desc1 - empPostburnin$desc2)
  empMeanCOMB <- mean(((empPostburnin$asc1 + empPostburnin$desc1) / 2) - ((empPostburnin$asc2 + empPostburnin$desc2) / 2))
  if(polyploidy == T){
    empMeanPOL <- empPostburnin$pol1 - empPostburnin$pol2
  }
  # calculate Delta R of each simulation
  Delta.R <- as.data.frame(matrix(data = NA, 
                                  nrow = nsim,
                                  ncol = 4))
  colnames(Delta.R) <- c("asc", "desc", "comb", "pol")
  # calculate empirical P
  for(i in 1:nsim){
    # get post burnin
    if(is.data.frame(simMCMC) == T){
      simPostburnin <- getPostBurnin(data = simMCMC,
                                     burn = burn)
    }else{
      simPostburnin <- getPostBurnin(data = simMCMC[[i]],
                                     burn = burn)
    }
    # isolate post burnin of each rate class
    pbrnASC1 <- simPostburnin$asc1
    pbrnASC2 <- simPostburnin$asc2
    pbrnDESC1 <- simPostburnin$desc1
    pbrnDESC2 <- simPostburnin$desc2
    pbrnCOMB1 <- (pbrnASC1 + pbrnDESC1) /2
    pbrnCOMB2 <- (pbrnASC2 + pbrnDESC2) /2
    # fill in DeltaR table
    Delta.R$asc[i] <- mean(abs(pbrnASC1 - pbrnASC2))
    Delta.R$desc[i] <- mean(abs(pbrnDESC1 - pbrnDESC2))
    Delta.R$comb[i] <- mean(abs(pbrnCOMB1 - pbrnCOMB2))
    if(polyploidy == T){
      pbrnPOL1 <- simPostburnin$pol1
      pbrnPOL2 <- simPostburnin$pol1
      Delta.R$pol[i] <- mean(abs(pbrnPOL1 - pbrnPOL2))
    }  
  }
  empericalPFission <- sum(Delta.R$asc > abs(empMeanASC)) / nsim
  empericalPFusion <- sum(Delta.R$desc > abs(empMeanDESC)) / nsim
  empericalPCombined <- sum(Delta.R$comb > abs(empMeanCOMB)) / nsim
  # print emperical p values
  print(paste("Emperical P value of Fusions:", empericalPFusion))
  print(paste("Emperical P value of Fissions:", empericalPFission))
  print(paste("Emperical P value of Combined aneuploidy:", empericalPCombined))
  if(polyploidy == T){
    empericalPPolyploidy <- sum(Delta.R$pol > abs(empMeanPOL))
    print(paste("Emperical P value of Polyploidy:", empericalPPolyploidy))
  }
  # get densities
  DensAsc <- density(Delta.R$asc)
  DensDesc <- density(Delta.R$desc)
  DensAnu <- density(Delta.R$comb)
  if(plot == T){
    par(mfcol = c(1,3))
    if(polyploidy == T){
      par(mfcol = c(1,4))
    }
    # fission
    plot(DensAsc, xlab = expression(paste("| ",Delta, "R"[Fission]," |")),
         main = "", cex.lab = 1.5, xlim  = c(0, max(c(DensAsc$x,abs(empMeanASC)))),
         lwd = 2)
    abline(v = mean(abs(empMeanASC)), col = "red", lwd = 2)
    mtext(text = "A", line = 0, outer = F, side = 3, adj = 0, cex = 1)
    # fusion
    plot(DensDesc, xlab = expression(paste("| ",Delta, "R"[Fusion]," |")),
         main = "", cex.lab = 1.5, 
         xlim  = c(0, max(c(DensDesc$x,abs(empMeanDESC)))), lwd = 2)
    abline(v = mean(abs(empMeanDESC)), col = "red", lwd = 2)
    mtext(text = "B", line = 0, outer = F, side = 3, adj = 0, cex = 1)
    # Combined aneuploidy
    plot(DensAnu, xlab = expression(paste("| ",Delta, "R"[Aneuploidy]," |")),
         main = "", cex.lab = 1.5,
         xlim  = c(0, max(c(DensAnu$x,abs(empMeanCOMB)))), lwd = 2)
    abline(v = mean(abs(empMeanCOMB)), col = "red", lwd = 2)
    mtext(text = "C", line = 0, outer = F, side = 3, adj = 0, cex = 1)
    if(polyploidy == T){
      DensPol <- density(Delta.R$pol)
      plot(DensPol, xlab = expression(paste("| ",Delta, "R"[polyploidy]," |")),
           main = "", cex.lab = 1.5, xlim  = c(0, max(c(DensPol$x,abs(empMeanPOL)))),
           lwd = 2)
      abline(v = mean(abs(empMeanPOL)), col = "red", lwd = 2)
      mtext(text = "D", line = 0, outer = F, side = 3, adj = 0, cex = 1)
    }
    par(mfcol = c(1,1))
  }
  # get p values
  Pvalues <-  list(EmpPvalueFusion = empericalPFusion,
                   EmpPvalueFission = empericalPFission,
                   EmpPvalueAneuploidy = empericalPCombined)
  if(polyploidy == T){
    Pvalues <-  list(EmpPvalueFusion = empericalPFusion,
                     EmpPvalueFission = empericalPFission,
                     EmpPvalueAneuploidy = empericalPCombined,
                     EmpPvaluePolyploidy = empericalPPolyploidy)
  }
  return(Pvalues)
} 

# This function will calculate HPD intervals ----
HPDcalc <- function(MCMC = NULL,
                    polyploidy = F,
                    nsim = 100,
                    plot = F,
                    burn = NULL){
  # make a table to get the HPD intervals for each MCMC we did
  HPD.fusions <- vector(mode = "numeric", length = nsim)
  HPD.fissions <- vector(mode = "numeric", length = nsim)
  HPD.combined <- vector(mode = "numeric", length = nsim)
  # fill in data
  for(i in 1:nsim){
    # get post burnin
    if(nsim == 1){
      empPostburnin <-  getPostBurnin(data = MCMC,burn = burn)
    }else{
      empPostburnin <-  getPostBurnin(data = MCMC[[i]],burn = burn)
    }
    # get data 
    empASC1 <- empPostburnin$asc1
    empASC2 <- empPostburnin$asc2
    empDESC1 <- empPostburnin$desc1
    empDESC2 <- empPostburnin$desc2
    empCOMB1 <- (empPostburnin$asc1 + empPostburnin$desc1) /2
    empCOMB2 <- (empPostburnin$asc2 + empPostburnin$desc2) /2
    # calculate HPD
    HPDasc <- HPDinterval(as.mcmc(empASC1 - empASC2))
    HPDdesc <- HPDinterval(as.mcmc(empDESC1 - empDESC2))
    HPDcomb <- HPDinterval(as.mcmc(empCOMB1 - empCOMB2))
    # calculate number of times HPD passes zero
    #fission
    if(HPDasc[1] < 0 & HPDasc[2] > 0){
      HPD.fissions[i] <- 0
    }else{
      HPD.fissions[i] <- 1
    }
    # fusion
    if(HPDdesc[1] < 0 & HPDdesc[2] > 0){
      HPD.fusions[i] <- 0
    }else{
      HPD.fusions[i] <- 1
    }
    # aneuploidy
    if(HPDcomb[1] < 0 & HPDcomb[2] > 0){
      HPD.combined[i] <- 0
    }else{
      HPD.combined[i] <- 1
    }
    # polyploidy
    if(polyploidy == T){
      HPD.polyploidy <- vector(mode = "numeric", length = nsim)
      #get data
      empPOL1 <- empPostburnin$asc1
      empPOL2 <- empPostburnin$asc2
      # calculate HPD
      HPDpol <- HPDinterval(as.mcmc(empASC1 - empASC2))
      # calculate number of times HPD passes zero
      if(HPDpol[1] < 0 & HPDpol[2] > 0){
        HPD.polyploidy[i] <- 0
      }else{
        HPD.polyploidy[i] <- 1
      }
    }
  }
  # get sum
  print(paste("Proportion of significant results for fissions:", sum(HPD.fissions)/ nsim))
  print(paste("Proportion of significant results for fusions:", sum(HPD.fusions)/ nsim))
  print(paste("Proportion of significant results for aneuploidy:", sum(HPD.combined)/ nsim))
  # get results
  results <- list(sigFission = sum(HPD.fissions)/ nsim,
                  sigFusion = sum(HPD.fusions)/ nsim,
                  sigAneuploidy =  sum(HPD.combined)/ nsim)
  return(results)
}

# This function will calculate empirical p value ----
# for a given tree and run chromeplus on that 
getEmpiricalP <- function(tree = NULL, 
                          chroms = NULL,
                          rng = NULL,
                          root.exp = NULL,
                          binary = NULL,
                          data = NULL,
                          nsim = NULL,
                          root = NULL,
                          args.lik = list(control = "ode",
                                          strict = F),
                          args.conlik = list(hyper = T,
                                             polyploidy = F,
                                             verbose = F,
                                             oneway = F,
                                             drop.poly=T,
                                             drop.demi=T,
                                             symmetric=F,
                                             nometa=F,
                                             meta="ARD"),
                          args.root = list(pi = "fitzjohn"),
                          args.MCMC = list(iter.temp = 20,
                                           iter = 100,
                                           prior = make.prior.exponential(r = .5),
                                           print.every=50),
                          burn = 0.5,
                          plot.lik = F,
                          plot.p = F){
  # run time complete script
  start_time_script <- as.numeric(Sys.time())
  #### checks ####
  #### checks ####
  #### ChromPlus ####
  # make a place holder to store ChromPlus output
  results <- vector(mode = "list",
                    length = nsim)
  # get Q matrix for binary state transitions
  Qmat <- matrix(data = 0, nrow = 2, ncol = 2)
  colnames(Qmat) <- rownames(Qmat) <- c(0,1)
  Qmat[1,2] <- mean(data$tran12)
  Qmat[2,1] <- mean(data$tran21)
  diag(Qmat) <- -rowSums(Qmat)
  # get the root probabilities of the binary trait
  root <- fitMk(tree = tree, 
                x = binary,
                fixedQ = Qmat,
                pi = args.root$pi)
  # simulate binary trait and run MCMC
  for(i in 1:nsim){
    # run time MCMC
    start_time_MCMC <- as.numeric(Sys.time())
    print(paste("simulation", i))
    check <- F
    while(check == F){
      # simulate binary traits
      sim.bin <- NULL
      sim.bin <- sim.character(tree = tree,
                               pars = Qmat,
                               model = "mkn",
                               x0 = as.numeric(sample(names(root$pi),
                                                      1,
                                                      prob = root$pi))+1)
      # change simulated binary characters back to zeros and ones
      sim.bin[sim.bin == 1] <- 0
      sim.bin[sim.bin == 2] <- 1
      if(sum(sim.bin) >= (length(sim.bin)*0.10) &
         sum(sim.bin) <= (length(sim.bin)*0.90)){
        check <- T
      }
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
    chrom.mat <- chromePlus::datatoMatrix(x = MCMC.dat,
                                          range = rng,
                                          hyper = args.conlik$hyper)
    # make the likelihood function
    lik <- NULL
    lik <- diversitree::make.mkn(tree = tree,
                                 states = chrom.mat,
                                 k = ncol(chrom.mat),
                                 strict = F,
                                 control = list(method="ode"))
    # constrain the likelihood function
    con.lik <- NULL
    con.lik <- chromePlus::constrainMkn(data = chrom.mat,
                                        lik = lik,
                                        hyper = args.conlik$hyper,
                                        polyploidy = args.conlik$polyploidy,
                                        verbose = args.conlik$verbose,
                                        oneway = args.conlik$oneway,
                                        constrain = list(drop.poly=args.conlik$drop.poly,
                                                         drop.demi=args.conlik$drop.demi,
                                                         symmetric = args.conlik$symmetric,
                                                         nometa = args.conlik$nometa,
                                                         meta = args.conlik$meta))
    # run the initial MCMC to get parameter values for w (only once)
    if(i == 1){
      temp <- NULL
      temp <- diversitree::mcmc(lik = con.lik,
                                x.init = runif(min=0, max=1,
                                               n=length(argnames(con.lik))),
                                prior = args.MCMC$prior,
                                nsteps = args.MCMC$iter.temp,
                                w = 1,
                                lower = rep(0,length(argnames(con.lik))),
                                print.every=5)
      # get values for w # TODO just do this first time through
      w <- NULL
      w <- diff(sapply(temp[(args.MCMC$iter.temp * 0.5 + 1):args.MCMC$iter.temp,
                            2:(length(argnames(con.lik))+1)],
                       quantile,
                       c(.05, .95)))
    }
    # run MCMC
    results[[i]] <-  diversitree::mcmc(lik = con.lik,
                                       x.init = runif(min=0, max=1,
                                                      n=length(argnames(con.lik))),
                                       nsteps = args.MCMC$iter,
                                       w = w,
                                       prior = args.MCMC$prior,
                                       lower = rep(0,length(argnames(con.lik))),
                                       print.every=args.MCMC$print.every)
    # print run time MCMC
    end_time_MCMC <- as.numeric(Sys.time())
    print(paste("run time MCMC =",  round(end_time_MCMC - start_time_MCMC,0), "seconds"))
  }
  if(nsim == 1){
    results <- results[[1]]
  }
  # plot likelihoods of the results
  if(plot.lik == T){
    plotlikMCMC(data = results, burn = burn)
  }
  #### End of ChromPlus ####
  #### calculation of empirical p value ####
  empP <- empiricalPcalc(empPostburnin = data,
                         simMCMC = results,
                         polyploidy = !args.conlik$drop.poly,
                         nsim = nsim,
                         plot = plot.p,
                         burn = burn)
  # run time complete script
  end_time_script <- as.numeric(Sys.time())
  print(paste("run time Total =",  round(end_time_script - start_time_script,0), "seconds"))
  # get results
  finalResults <- list(MCMC = results,
                       empiricalP = empP)
  return(finalResults)
}
# This function will calculate empirical p value in multicore ----
# for a given tree and run chromeplus on that 
getEmpiricalPMC <- function(tree = NULL, 
                            chroms = NULL,
                            rng = NULL,
                            root.exp = NULL,
                            binary = NULL,
                            data = NULL,
                            nsim = NULL,
                            root = NULL,
                            args.lik = list(control = "ode",
                                            strict = F),
                            args.conlik = list(hyper = T,
                                               polyploidy = F,
                                               verbose = F,
                                               oneway = F,
                                               drop.poly=T,
                                               drop.demi=T,
                                               symmetric=F,
                                               nometa=F,
                                               meta="ARD"),
                            args.root = list(pi = "fitzjohn"),
                            args.MCMC = list(iter.temp = 20,
                                             iter = 100,
                                             prior = make.prior.exponential(r = .5),
                                             print.every=50),
                            burn = 0.5,
                            plot.lik = F,
                            plot.p = F,
                            nclust = 2){
  # run time complete script
  start_time_script <- as.numeric(Sys.time())
  #### checks ####
  #### checks ####
  #### ChromPlus ####
  # define number of clusters for parallel computing
  NumberOfClusters <- nclust
  cl <- makeCluster(NumberOfClusters, outfile = "")
  registerDoSNOW(cl)
  # get Q matrix for binary state transitions
  Qmat <- matrix(data = 0, nrow = 2, ncol = 2)
  colnames(Qmat) <- rownames(Qmat) <- c(0,1)
  Qmat[1,2] <- mean(data$tran12)
  Qmat[2,1] <- mean(data$tran21)
  diag(Qmat) <- -rowSums(Qmat)
  # get the root probabilities of the binary trait
  root <- fitMk(tree = tree, 
                x = binary,
                fixedQ = Qmat,
                pi = args.root$pi)
  # simulate binary trait and run MCMC
  results <- foreach(i = 1:nsim, .verbose = T, .packages = c("ape","diversitree", "chromePlus","phytools","maps")) %dopar% {
    # run time MCMC
    start_time_MCMC <- as.numeric(Sys.time())
    print(paste("simulation", i))
    check <- F
    while(check == F){
      # simulate binary traits
      sim.bin <- NULL
      sim.bin <- sim.character(tree = tree,
                               pars = Qmat,
                               model = "mkn",
                               x0 = as.numeric(sample(names(root$pi),
                                                      1,
                                                      prob = root$pi))+1)
      # change simulated binary characters back to zeros and ones
      sim.bin[sim.bin == 1] <- 0
      sim.bin[sim.bin == 2] <- 1
      if(sum(sim.bin) >= (length(sim.bin)*0.10) &
         sum(sim.bin) <= (length(sim.bin)*0.90)){
        check <- T
      }
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
    chrom.mat <- chromePlus::datatoMatrix(x = MCMC.dat,
                                          range = rng,
                                          hyper = args.conlik$hyper)
    # make the likelihood function
    lik <- NULL
    lik <- diversitree::make.mkn(tree = tree,
                                 states = chrom.mat,
                                 k = ncol(chrom.mat),
                                 strict = F,
                                 control = list(method="ode"))
    # constrain the likelihood function
    con.lik <- NULL
    con.lik <- chromePlus::constrainMkn(data = chrom.mat,
                                        lik = lik,
                                        hyper = args.conlik$hyper,
                                        polyploidy = args.conlik$polyploidy,
                                        verbose = args.conlik$verbose,
                                        oneway = args.conlik$oneway,
                                        constrain = list(drop.poly=args.conlik$drop.poly,
                                                         drop.demi=args.conlik$drop.demi,
                                                         symmetric = args.conlik$symmetric,
                                                         nometa = args.conlik$nometa,
                                                         meta = args.conlik$meta))
    # run the initial MCMC to get parameter values for w (only once)
    temp <- NULL
    temp <- diversitree::mcmc(lik = con.lik,
                              x.init = runif(min=0, max=1,
                                             n=length(argnames(con.lik))),
                              prior = args.MCMC$prior,
                              nsteps = args.MCMC$iter.temp,
                              w = 1,
                              lower = rep(0,length(argnames(con.lik))),
                              print.every=5)
    # get values for w # TODO just do this first time through
    w <- NULL
    w <- diff(sapply(temp[(args.MCMC$iter.temp * 0.5 + 1):args.MCMC$iter.temp,
                          2:(length(argnames(con.lik))+1)],
                     quantile,
                     c(.05, .95)))
    # run MCMC
    results <-  diversitree::mcmc(lik = con.lik,
                                  x.init = runif(min=0, max=1,
                                                 n=length(argnames(con.lik))),
                                  nsteps = args.MCMC$iter,
                                  w = w,
                                  prior = args.MCMC$prior,
                                  lower = rep(0,length(argnames(con.lik))),
                                  print.every=args.MCMC$print.every)
    # print run time MCMC
    end_time_MCMC <- as.numeric(Sys.time())
    print(paste("run time MCMC =",  round(end_time_MCMC - start_time_MCMC,0), "seconds"))
    
    final_result <- results
  }
  # stop the cluster
  stopCluster(cl)
  
  if(nsim == 1){
    results <- results[[1]]
  }
  # plot likelihoods of the results
  if(plot.lik == T){
    plotlikMCMC(data = results, burn = burn)
  }
  #### End of ChromPlus ####
  #### calculation of empirical p value ####
  empP <- empiricalPcalc(empPostburnin = data,
                         simMCMC = results,
                         polyploidy = !args.conlik$drop.poly,
                         nsim = nsim,
                         plot = plot.p,
                         burn = burn)
  # run time complete script
  end_time_script <- as.numeric(Sys.time())
  print(paste("run time Total =",  round(end_time_script - start_time_script,0), "seconds"))
  # get results
  finalResults <- list(MCMC = results,
                       empiricalP = empP)
  return(finalResults)
}