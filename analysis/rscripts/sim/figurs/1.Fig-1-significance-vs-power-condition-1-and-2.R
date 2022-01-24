# Terrence Sylvester
# pradakshanas@gmail.com
# 10 January 2022

# load libraries
library(coda)
library(stringr)

#### calculate significance ####
# get the list of files
files <- dir("../results/post.burnin/")
# only keep files that we want
files <- files[grep("post.burnin", files)]
files <- files[grep(paste("condition.1", sep = ""), files)]
# container to store p values
Sigtab <- as.data.frame(matrix(data = NA,
                               nrow = 15,
                               ncol = 6))
colnames(Sigtab) <- c("condition", "ntip","rr", "FissionSignificance", "FusionSignificance", "MeanRateSignificance")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# read in data
for(i in 1:length(files)){
  dat <- read.csv(paste("../results/post.burnin/", files[i], sep = ""))
  # make tables to hold rates in each simulation
  tab <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 3))
  # label columns
  colnames(tab) <- c("FissionPassZero", "FusionPassZero", "MeanRatepassZero")
  for(j in 1:100){
    # fission
    # calculate HPD interval
    HPD.fis <- HPDinterval(as.mcmc(dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$asc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fis[1] < 0 & HPD.fis[2] > 0){
      tab$FissionPassZero[j] <- 0
    }else{
      tab$FissionPassZero[j] <- 1
    }
    # fusion
    # calculate HPD interval
    HPD.fus <- HPDinterval(as.mcmc(dat$desc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$desc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fus[1] < 0 & HPD.fus[2] > 0){
      tab$FusionPassZero[j] <- 0
    }else{
      tab$FusionPassZero[j] <- 1
    }
    # combined
    # calculate HPD interval
    HPD.comb <- HPDinterval(as.mcmc((dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc2[pb.seq[j]:(pb.seq[j] + 49)])/2  - (dat$asc1[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc1[pb.seq[j]:(pb.seq[j] + 49)])/2))
    if(HPD.comb[1] < 0 & HPD.comb[2] > 0){
      tab$MeanRatepassZero[j] <- 0
    }else{
      tab$MeanRatepassZero[j] <- 1
    }
  }
  
  Sigtab$condition[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                                "tree50." = "",
                                                                "tree100." = "",
                                                                "tree150." = "",
                                                                "tree200." = "",
                                                                "tree250." = "",
                                                                "condition." = "",
                                                                ".rr.2.csv" = "",
                                                                ".rr.5.csv" = "",
                                                                ".rr.10.csv" = "")))
  
  Sigtab$ntip[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                           "condition.1." = "",
                                                           "condition.2." = "",
                                                           "condition.3." = "",
                                                           "condition.4." = "",
                                                           "condition.5." = "",
                                                           "condition.6." = "",
                                                           "tree" = "",
                                                           ".rr.2.csv" = "",
                                                           ".rr.5.csv" = "",
                                                           ".rr.10.csv" = "")))
  
  Sigtab$rr[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                         "condition.1." = "",
                                                         "condition.2." = "",
                                                         "condition.3." = "",
                                                         "condition.4." = "",
                                                         "condition.5." = "",
                                                         "condition.6." = "",
                                                         "tree50." = "",
                                                         "tree100." = "",
                                                         "tree150." = "",
                                                         "tree200." = "",
                                                         "tree250." = "",
                                                         "rr." = "",
                                                         ".csv" = "")))
  
  Sigtab$FissionSignificance[i] <- sum(tab$FissionPassZero) / 100
  Sigtab$FusionSignificance[i] <- sum(tab$FusionPassZero) / 100
  Sigtab$MeanRateSignificance[i] <- sum(tab$MeanRatepassZero) / 100
}
# re order ptab by tip number
Sigtab <- Sigtab[order(Sigtab$ntip),]
# keep onlu SigTab
rm(list = ls()[-c(which(ls() == "Sigtab"))]) 

#### calculate power ####
# get the list of files
files <- dir("../results/post.burnin/")
# only keep files that we want
files <- files[grep("post.burnin", files)]
files <- files[grep(paste("condition.2", sep = ""), files)]
# container to store p values
Ptab <- as.data.frame(matrix(data = NA,
                             nrow = 15,
                             ncol = 6))
colnames(Ptab) <- c("condition", "ntip","rr", "FissionSignificance", "FusionSignificance", "MeanRateSignificance")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# read in data
for(i in 1:length(files)){
  dat <- read.csv(paste("../results/post.burnin/", files[i], sep = ""))
  # make tables to hold rates in each simulation
  tab <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 3))
  # label columns
  colnames(tab) <- c("FissionPassZero", "FusionPassZero", "MeanRatepassZero")
  for(j in 1:100){
    # fission
    # calculate HPD interval
    HPD.fis <- HPDinterval(as.mcmc(dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$asc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fis[1] < 0 & HPD.fis[2] > 0){
      tab$FissionPassZero[j] <- 0
    }else{
      tab$FissionPassZero[j] <- 1
    }
    # fusion
    # calculate HPD interval
    HPD.fus <- HPDinterval(as.mcmc(dat$desc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$desc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fus[1] < 0 & HPD.fus[2] > 0){
      tab$FusionPassZero[j] <- 0
    }else{
      tab$FusionPassZero[j] <- 1
    }
    # combined
    # calculate HPD interval
    HPD.comb <- HPDinterval(as.mcmc((dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc2[pb.seq[j]:(pb.seq[j] + 49)])/2  - (dat$asc1[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc1[pb.seq[j]:(pb.seq[j] + 49)])/2))
    if(HPD.comb[1] < 0 & HPD.comb[2] > 0){
      tab$MeanRatepassZero[j] <- 0
    }else{
      tab$MeanRatepassZero[j] <- 1
    }
  }
  
  Ptab$condition[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                              "tree50." = "",
                                                              "tree100." = "",
                                                              "tree150." = "",
                                                              "tree200." = "",
                                                              "tree250." = "",
                                                              "condition." = "",
                                                              ".rr.2.csv" = "",
                                                              ".rr.5.csv" = "",
                                                              ".rr.10.csv" = "")))
  
  Ptab$ntip[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                         "condition.1." = "",
                                                         "condition.2." = "",
                                                         "condition.3." = "",
                                                         "condition.4." = "",
                                                         "condition.5." = "",
                                                         "condition.6." = "",
                                                         "tree" = "",
                                                         ".rr.2.csv" = "",
                                                         ".rr.5.csv" = "",
                                                         ".rr.10.csv" = "")))
  
  Ptab$rr[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                       "condition.1." = "",
                                                       "condition.2." = "",
                                                       "condition.3." = "",
                                                       "condition.4." = "",
                                                       "condition.5." = "",
                                                       "condition.6." = "",
                                                       "tree50." = "",
                                                       "tree100." = "",
                                                       "tree150." = "",
                                                       "tree200." = "",
                                                       "tree250." = "",
                                                       "rr." = "",
                                                       ".csv" = "")))
  
  Ptab$FissionSignificance[i] <- sum(tab$FissionPassZero) / 100
  Ptab$FusionSignificance[i] <- sum(tab$FusionPassZero) / 100
  Ptab$MeanRateSignificance[i] <- sum(tab$MeanRatepassZero) / 100
}
# re order ptab by tip number
Ptab <- Ptab[order(Ptab$ntip),]
# keep onlu SigTab and Ptab
rm(list = ls()[-c(which(ls() %in% c("Sigtab", "Ptab")))]) 
# plot
{
  par(mfrow = c(2,3))
  #set colours and point types
  col <- viridis::viridis(3, alpha = 0.5, end = 0.8)
  
  ## Type I error rate Fission
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,0.1),
       ylab = "Type I error rate",
       xlab = "Number of tips",
       main = "Fission") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          y = Sigtab$FissionSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
  # add 0.05 p value cut off
  abline(h = 0.05, col = "red", lty = 2)
  
  ## Type I error rate Fusion
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,0.1),
       ylab = "Type I error rate",
       xlab = "Number of tips",
       main = "Fusion") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          y = Sigtab$FusionSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
  # add 0.05 p value cut off
  abline(h = 0.05, col = "red", lty = 2)
  
  ## Type I error rate MeanRateCombined
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,0.1),
       ylab = "Type I error rate",
       xlab = "Number of tips",
       main = "Combination") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Sigtab$ntip[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          y = Sigtab$MeanRateSignificance[Sigtab$rr == sort(unique(Sigtab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
  # add legend
  legend("topleft",inset=.02,
         legend=c("2",
                  "5",
                  "10"),
         col=col, 
         lty=1,
         cex=1,
         lwd = 2,title = "Rates ratio")
  
  # add 0.05 p value cut off
  abline(h = 0.05, col = "red", lty = 2)
  
  ## Type II error rate fission
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,1),
       ylab = "Type II error rate",
       xlab = "Number of tips") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
          y = Ptab$FissionSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
  
  ## Type II error rate fusion
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,1),
       ylab = "Type II error rate",
       xlab = "Number of tips") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
          y = Ptab$FusionSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
  
  ## Type II error rate MeanRateCombined
  plot(x = NULL,
       y = NULL,
       xlim = c(50,250),
       ylim = c(0,1),
       ylab = "Type II error rate",
       xlab = "Number of tips") 
  
  for(i in 1:3){
    # plot condition 2
    lines(x = Ptab$ntip[Ptab$rr == sort(unique(Ptab$rr))[i]],
          y = Ptab$MeanRateSignificance[Ptab$rr == sort(unique(Ptab$rr))[i]],
          pch = 16,
          type = "o", lwd = 2,
          col = col[i])
    
  }
}
