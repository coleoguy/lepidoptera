# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# load libraries
library(coda)
# get helper functions
source("../rscripts/functions.R")
# list MCMC results
files <- dir("../results/simData-MCMC/")
files <- files[grep("cond", files)]
# make a data frame to hold final results
dat <- as.data.frame(matrix(data = NA,
                            nrow = length(files),
                            ncol = 4))
colnames(dat) <- c("sim", "sigFission", "sigFusion", "sigAneuploidy")
# run loop to get power at each condition
for(i in 1:length(files)){
  # load MCMC results
  load(paste("../results/simData-MCMC/",files[i], sep = ""))
  # plot likelihood
  plotlikMCMC(data = results, burn = 0.5)
  # get number of significant results
  sigResults <- HPDcalc(MCMC = results,
                        polyploidy = F,
                        nsim = 100,
                        plot = F,
                        burn = 0.5)
  # fill in the table
  dat$sim[i] <- files[i]
  dat$sigFission[i] <- sigResults$sigFission
  dat$sigFusion[i] <- sigResults$sigFusion
  dat$sigAneuploidy[i] <- sigResults$sigAneuploidy
  # remove unwanted results
  rm(results, pbrn, sigResults)
}
# make new columns to get ntip, rr and condition for plotting
dat$tips <- rep(rep(c(100,200,50), each = 4), 3)
dat$rr <- rep(c(1,10,2,5),9)
dat$cond <- rep(c(1,2,3), each = 12)
# sort data table by tip number
dat <- dat[order(dat$tips),]
# plot
par(mfrow = c(3,3))
for(i in 1:3){
  # Fission
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigFission[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigFission[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigFission[dat$cond == i & dat$rr == 5],type = "o",col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigFission[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = paste("Condition", i),side = 3,adj = 0,
        cex = 0.7)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 0.7)
  # if(i == 1){
  # # # add legend
  # legend("topleft",inset=.02,
  #        legend=c("1",
  #                 "2",
  #                 "5",
  #                 "10"),
  #        col=c("gray","#e66101" ,"#fdb863" ,"#5e3c99"),
  #        lty=1,
  #        cex=1,
  #        lwd = 2,title = "Rates ratio")
  # }
  # Fusuin
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigFusion[dat$cond == i & dat$rr == 1],type = "o", col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigFusion[dat$cond == i & dat$rr == 2],type = "o", col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigFusion[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigFusion[dat$cond == i & dat$rr == 10],type = "o", col = "#5e3c99", pch = 16)
  mtext(text = "Fusion",side = 3,adj = 1,
        cex = 0.7)
  # Aneuploidy
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1))
  lines(x = dat$tips[dat$cond == i & dat$rr == 1], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 1],type = "o",col = "gray", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 2], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 2],type = "o",col = "#e66101", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 5], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 5],type = "o", col = "#fdb863", pch = 16)
  lines(x = dat$tips[dat$cond == i & dat$rr == 10], y = dat$sigAneuploidy[dat$cond == i & dat$rr == 10],type = "o",col = "#5e3c99", pch = 16)
  mtext(text = "Aneuploidy",side = 3,adj = 1,
        cex = 0.7)
  
}


