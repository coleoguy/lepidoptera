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
# save results
write.csv(dat,"../results/power.old.method.csv", row.names = F)
