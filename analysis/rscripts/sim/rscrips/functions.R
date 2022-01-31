# This function will perform the MCMC ----
runMCMC <- function(tree = NULL,
                    chroms = NULL,
                    binary = NULL,
                    hyper = T,
                    polyploidy = F,
                    verbose = F,
                    oneway = F,
                    drop.poly=T,
                    drop.demi=T,
                    iter.temp = 20,
                    iter = 100,
                    prior = make.prior.exponential(r = .5),
                    print.every=50){
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
  chrom.mat <- NULL
  chrom.mat <- datatoMatrix(x = MCMC.dat,
                            range = rng,
                            hyper = hyper)
  
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
                          hyper = hyper,
                          polyploidy = polyploidy,
                          verbose = verbose,
                          oneway = oneway,
                          constrain = list(drop.poly=drop.poly,
                                           drop.demi=drop.demi))
  
  # run the initial MCMC to get parameter values for w
  temp <- NULL
  temp <- diversitree::mcmc(lik = con.lik,
                            x.init = runif(min=0, max=1,
                                           n=length(argnames(con.lik))),
                            prior = prior,
                            nsteps = iter.temp,
                            w = 1,
                            lower = rep(0,length(argnames(con.lik))),
                            print.every = 5)
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
                               lower = rep(0,length(argnames(con.lik))),
                               print.every = print.every)
  end_time <- as.numeric(Sys.time())
  print(paste("run time =",  round(end_time - start_time,0), "seconds"))
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
  }
  
  if(is.list(data) == T){
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
  if(is.list(data) == T){
    y <- c()
    y <- data[[1]][((nrow(data[[1]])*burn)+1):nrow(data[[1]]),]
    for(i in 2:length(data)){
      y <- rbind(y,data[[i]][((nrow(data[[1]])*burn)+1):nrow(data[[1]]),])
    }
    pbrn <- y
  }
  if(is.data.frame(data) == T){
    pbrn <- data[[((nrow(data)*burn)+1):nrow(data),]]
  }
  return(pbrn)
}

# This function will make transision matrix  ----
# get the transistion matrix
# tmat <- matrix(data = 0, nrow = 2, ncol = 2)
# colnames(tmat) <- rownames(tmat) <- c(0,1)
# # fill tmat
# tmat[1,2] <- mean(postBurnResults$tran12)
# tmat[2,1] <- mean(postBurnResults$tran21)
# diag(tmat) <- -rowSums(tmat)



# This function will simulate a set of binary characters for a given tree and run MCMC on that ----
simBinMCMC <- function(tree = NULL,
                       chroms = NULL,
                       binary = NULL,
                       tmat = NULL,
                       nsim = NULL,
                       hyper = T,
                       polyploidy = F,
                       verbose = F,
                       oneway = F,
                       drop.poly=T,
                       drop.demi=T,
                       iter.temp = 20,
                       iter = 100,
                       prior = make.prior.exponential(r = .5),
                       print.every=50){
  
  results <- vector(mode = "list", length = nsim)
  # # get the root probabilities of the binary trait
  root <- fitMk(tree = tree, 
                x = binary,
                fixedQ = tmat,
                pi = "fitzjohn")
  # run time
  start_time <- as.numeric(Sys.time())
  for(i in 1:nsim){
    print(paste("iteration", i))
    
    check <- F
    while(check == F){
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
                            hyper = hyper,
                            polyploidy = polyploidy,
                            verbose = verbose,
                            oneway = oneway,
                            constrain = list(drop.poly=drop.poly,
                                             drop.demi=drop.demi))
    
    # run the initial MCMC to get parameter values for w
    temp <- NULL
    temp <- diversitree::mcmc(lik = con.lik,
                              x.init = runif(min=0, max=1,
                                             n=length(argnames(con.lik))),
                              prior = prior,
                              nsteps = iter.temp,
                              w = 1,
                              lower = rep(0,length(argnames(con.lik))),
                              print.every=5)
    # get values for w
    w <- NULL
    w <- diff(sapply(temp[11:20, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
    # run MCMC
    results[[i]] <-  diversitree::mcmc(lik = con.lik,
                                       x.init = runif(min=0, max=1,
                                                      n=length(argnames(con.lik))),
                                       nsteps = iter,
                                       w = w,
                                       prior = prior,
                                       lower = rep(0,length(argnames(con.lik))),
                                       print.every=print.every)
    
  }
  end_time <- as.numeric(Sys.time())
  print(paste("run time =",  round(end_time - start_time,0), "seconds"))
  return(results)
}

# This function will calculate empirical p value ----
empiricalPcalc <- function(empPostburnin = NULL,
                           simPostburnin = NULL,
                           polyploidy = F,
                           nsim = 100,
                           plot = F){
  
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
  # this will be used to get the post burnin from each MCMC run
  pb.seq <- seq(from = 1, to = 5000, by = 50)
  for(i in 1:nsim){
    pbrnASC1 <- simPostburnin$asc1[pb.seq[i]:(pb.seq[i]+49)]
    pbrnASC2 <- simPostburnin$asc2[pb.seq[i]:(pb.seq[i]+49)]
    pbrnDESC1 <- simPostburnin$desc1[pb.seq[i]:(pb.seq[i]+49)]
    pbrnDESC2 <- simPostburnin$desc2[pb.seq[i]:(pb.seq[i]+49)]
    pbrnCOMB1 <- (pbrnASC1 + pbrnDESC1) /2
    pbrnCOMB2 <- (pbrnASC2 + pbrnDESC2) /2
    
    Delta.R$asc[i] <- mean(abs(pbrnASC1 - pbrnASC2))
    Delta.R$desc[i] <- mean(abs(pbrnDESC1 - pbrnDESC2))
    Delta.R$comb[i] <- mean(abs(pbrnCOMB1 - pbrnCOMB2))
    
    
    if(polyploidy == T){
      pbrnPOL1 <- simPostburnin$pol1[pb.seq[i]:(pb.seq[i]+49)]
      pbrnPOL2 <- simPostburnin$pol1[pb.seq[i]:(pb.seq[i]+49)]
      Delta.R$pol[i] <- mean(abs(pbrnPOL1 - pbrnPOL2))
    }  
  }
  empericalPFission <- sum(Delta.R$asc > abs(empMeanASC)) / nsim
  empericalPFusion <- sum(Delta.R$desc > abs(empMeanDESC)) / nsim
  empericalPCombined <- sum(Delta.R$comb > abs(empMeanCOMB)) / nsim
  
  print(paste("Emperical P value of Fusions:", empericalPFusion))
  print(paste("Emperical P value of Fissions:", empericalPFission))
  print(paste("Emperical P value of Combined aneuploidy:", empericalPCombined))
  
  if(polyploidy == T){
    empericalPPolyploidy <- sum(Delta.R$pol > abs(empMeanPOL))
    print(paste("Emperical P value of Polyploidy:", empericalPPolyploidy))
  }
  
  if(plot == T){
    par(mfcol = c(1,3))
    if(polyploidy == T){
      par(mfcol = c(1,4))
    }
    # fission
    hist(Delta.R$asc,
         xlab = expression(paste("| ",Delta, "R"[Fission]," |")),
         main = "",
         breaks = 25,
         cex.lab = 1.5,
         lwd = 2,
         xlim  = c(0, max(c(Delta.R$asc,abs(empMeanASC)))))
    abline(v = mean(abs(empMeanASC)),
           col = "red",
           lwd = 2)
    mtext(text = "A",
          line = 0,
          outer = F,
          side = 3,
          adj = 0,
          cex = 1)
    
    # fusion
    hist(Delta.R$desc,
         xlab = expression(paste("| ",Delta, "R"[Fusion]," |")),
         main = "",
         breaks = 25,
         cex.lab = 1.5,
         lwd = 2,
         xlim  = c(0, max(c(Delta.R$desc,abs(empMeanDESC)))))
    abline(v = mean(abs(empMeanDESC)),
           col = "red",
           lwd = 2)
    mtext(text = "B",
          line = 0,
          outer = F,
          side = 3,
          adj = 0,
          cex = 1)
    
    # Combined aneuploidy
    hist(Delta.R$comb,
         xlab = expression(paste("| ",Delta, "R"[Aneuploidy]," |")),
         main = "",
         breaks = 25,
         cex.lab = 1.5,
         lwd = 2,
         xlim  = c(0, max(c(Delta.R$comb,abs(empMeanCOMB)))))
    abline(v = mean(abs(empMeanCOMB)),
           col = "red",
           lwd = 2)
    mtext(text = "C",
          line = 0,
          outer = F,
          side = 3,
          adj = 0,
          cex = 1)
    
    if(polyploidy == T){
      # Combined Polyploidy
      hist(Delta.R$pol,
           xlab = expression(paste("| ",Delta, "R"[polyploidy]," |")),
           main = "",
           breaks = 25,
           cex.lab = 1.5,
           lwd = 2)
      abline(v = mean(abs(empMeanPOL)),
             col = "red",
             lwd = 2)
      mtext(text = "D",
            line = 0,
            outer = F,
            side = 3,
            adj = 0,
            cex = 1)
    }
  }
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
HPDcalc <- function(empPostburnin = NULL,
                    polyploidy = F,
                    nsim = 100,
                    plot = F){
  # make a table to get the HPD intervals for each MCMC we did
  HPD.fusions <- vector(mode = "numeric", length = nsim)
  HPD.fissions <- vector(mode = "numeric", length = nsim)
  HPD.combined <- vector(mode = "numeric", length = nsim)
  
  # this will be used to get the post burnin from each MCMC run
  pb.seq <- seq(from = 1, to = 5000, by = 50)
  
  # fill in data
  for(i in 1:nsim){
    
    # get data 
    empASC1 <- empPostburnin$asc1[pb.seq[i]:(pb.seq[i]+49)]
    empASC2 <- empPostburnin$asc2[pb.seq[i]:(pb.seq[i]+49)]
    empDESC1 <- empPostburnin$desc1[pb.seq[i]:(pb.seq[i]+49)]
    empDESC2 <- empPostburnin$desc2[pb.seq[i]:(pb.seq[i]+49)]
    empCOMB1 <- (empPostburnin$asc1[pb.seq[i]:(pb.seq[i]+49)] + empPostburnin$desc1[pb.seq[i]:(pb.seq[i]+49)]) /2
    empCOMB2 <- (empPostburnin$asc2[pb.seq[i]:(pb.seq[i]+49)] + empPostburnin$desc2[pb.seq[i]:(pb.seq[i]+49)]) /2
    
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
      empPOL1 <- empPostburnin$asc1[pb.seq[i]:(pb.seq[i]+49)]
      empPOL2 <- empPostburnin$asc2[pb.seq[i]:(pb.seq[i]+49)]
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
  
  results <- list(sigFission = sum(HPD.fissions)/ nsim,
                  sigFusion = sum(HPD.fusions)/ nsim,
                  sigAneuploidy =  sum(HPD.combined)/ nsim)
  return(results)
  
  # sum(HPD.tab.fissions$`pass-zero`)
  
  # # plot data
  # par(mfrow = c(2,3))
  # #set colours and point types
  # col <- viridis::viridis(3, alpha = 0.5, end = 0.8)
  # 
  # ## Type I error rate Fission
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,0.1),
  #      ylab = "Type I error rate",
  #      xlab = "Number of tips",
  #      main = "Fission") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         y = Sigtab$FissionSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
  # }
  # # add 0.05 p value cut off
  # abline(h = 0.05, col = "red", lty = 2)
  # 
  # ## Type I error rate Fusion
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,0.1),
  #      ylab = "Type I error rate",
  #      xlab = "Number of tips",
  #      main = "Fusion") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         y = Sigtab$FusionSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
  # }
  # # add 0.05 p value cut off
  # abline(h = 0.05, col = "red", lty = 2)
  # 
  # ## Type I error rate MeanRateCombined
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,0.1),
  #      ylab = "Type I error rate",
  #      xlab = "Number of tips",
  #      main = "Combination") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         y = Sigtab$MeanRateSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
  # }
  # # add legend
  # legend("topleft",inset=.02,
  #        legend=c("2",
  #                 "5",
  #                 "10"),
  #        col=col, 
  #        lty=1,
  #        cex=1,
  #        lwd = 2,title = "Rates ratio")
  # 
  # # add 0.05 p value cut off
  # abline(h = 0.05, col = "red", lty = 2)
  # 
  # ## Type II error rate fission
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,1),
  #      ylab = "Type II error rate",
  #      xlab = "Number of tips") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         y = Ptab$FissionSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
  # }
  # 
  # ## Type II error rate fusion
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,1),
  #      ylab = "Type II error rate",
  #      xlab = "Number of tips") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         y = Ptab$FusionSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
  # }
  # 
  # ## Type II error rate MeanRateCombined
  # plot(x = NULL,
  #      y = NULL,
  #      xlim = c(50,250),
  #      ylim = c(0,1),
  #      ylab = "Type II error rate",
  #      xlab = "Number of tips") 
  # 
  # for(i in 1:3){
  #   # plot condition 2
  #   lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         y = Ptab$MeanRateSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
  #         pch = 16,
  #         type = "o", lwd = 2,
  #         col = col[i])
  #   
}

