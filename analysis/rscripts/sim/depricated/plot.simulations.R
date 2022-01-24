# Terrence Sylvester
# pradakshanas@gmail.com
# 10 January 2022

# load libraries
library(coda)
library(stringr)

#### conditions ####
# condition 1 : No signal for binary trait and no tip violation
# condition 2 : Signal for binary trait and no tip violation
# condition 3 : No signal for binary trait and Single tip violation
# condition 4 : Signal for binary trait and Single tip violation
# condition 5 : No signal for binary trait and tip violators increase up to 10%
# condition 6 : Signal for binary trait and tip violators increase up to 10%

#### simulation parameters ####

# simulation parameters for condition 1,3, and 5
# rate of fusion at state 1       > 0.75
# rate of fission at state 1      > 0.75
# rate of fusion at state 0       > 0.75
# rate of fission at state 0      > 0.75
# rate of polyploidy at state 1   > 0
# rate of demiploidy at state 1   > 0
# rate of polyploidy at state 0   > 0
# rate of demiploidy at state 0   > 0
# transision from state 0 to 1    > 0.5
# transision from state 1 to 0    > 0.5
# root chromosome number          > 10
# root state                      > 0

# simulation parameters for condition 2,4, and 6
# rate of fusion at state 1       > 0.75 * rates ratio
# rate of fission at state 1      > 0.75 * rates ratio
# rate of fusion at state 0       > 0.75
# rate of fission at state 0      > 0.75
# rate of polyploidy at state 1   > 0
# rate of demiploidy at state 1   > 0
# rate of polyploidy at state 0   > 0
# rate of demiploidy at state 0   > 0
# transision from state 0 to 1    > 0.5
# transision from state 1 to 0    > 0.5
# root chromosome number          > 10
# root state                      > 0

#### process data ####
# get the list of files
files <- dir()
# only keep files that we want
files <- files[grep("post.burnin", files)]
# get the p value
ptab <- as.data.frame(matrix(data = NA,
                             nrow = 90,
                             ncol = 9))
colnames(ptab) <- c("condition",
                    "ntip",
                    "rr",
                    "pvalue.fis",
                    "pvalue.fus",
                    "pvalue.comb",
                    "passZeroFis",
                    "passZeroFus",
                    "passZeroComb")
# this will be used to get the post burnin from each MCMC run
pb.seq <- seq(from = 1, to = 5000, by = 50)
# read in data
for(i in 1:length(files)){
  dat <- read.csv(files[i])
  # make tables to hold rates in each simulation
  tab <- as.data.frame(matrix(data = NA, nrow = 100, ncol = 6))
  # label columns
  colnames(tab) <- c("D.fis", "D.fis.zero","D.fus", "D.fus.zero", "D.comb", "D.comb.zero")
  for(j in 1:100){
    # fission
    tab$D.fis[j] <- mean(dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$asc1[pb.seq[j]:(pb.seq[j] + 49)])
    HPD.fis <- HPDinterval(as.mcmc(dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$asc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fis[1] < 0 & HPD.fis[2] > 0){
      tab$D.fis.zero[j] <- 0
    }else{
      tab$D.fis.zero[j] <- 1
    }
    # fusion
    tab$D.fus[j] <- mean(dat$desc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$desc1[pb.seq[j]:(pb.seq[j] + 49)])
    HPD.fus <- HPDinterval(as.mcmc(dat$desc2[pb.seq[j]:(pb.seq[j] + 49)] - dat$desc1[pb.seq[j]:(pb.seq[j] + 49)]))
    if(HPD.fus[1] < 0 & HPD.fus[2] > 0){
      tab$D.fus.zero[j] <- 0
    }else{
      tab$D.fus.zero[j] <- 1
    }
    # combined
    tab$D.comb[j] <- mean((dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc2[pb.seq[j]:(pb.seq[j] + 49)])/2  - (dat$asc1[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc1[pb.seq[j]:(pb.seq[j] + 49)])/2)
    HPD.comb <- HPDinterval(as.mcmc((dat$asc2[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc2[pb.seq[j]:(pb.seq[j] + 49)])/2  - (dat$asc1[pb.seq[j]:(pb.seq[j] + 49)] + dat$desc1[pb.seq[j]:(pb.seq[j] + 49)])/2))
    if(HPD.comb[1] < 0 & HPD.comb[2] > 0){
      tab$D.comb.zero[j] <- 0
    }else{
      tab$D.comb.zero[j] <- 1
    }
  }
  ptab$condition[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
                                                              "tree50." = "",
                                                              "tree100." = "",
                                                              "tree150." = "",
                                                              "tree200." = "",
                                                              "tree250." = "",
                                                              "condition." = "",
                                                              ".rr.2.csv" = "",
                                                              ".rr.5.csv" = "",
                                                              ".rr.10.csv" = "")))
  ptab$ntip[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
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
  ptab$rr[i] <- as.numeric(str_replace_all(files[i], c("post.burnin." = "",
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
  if(ptab$condition[i] %in% c(2,4,6) & ptab$rr[i] == 2){
    ptab$pvalue.fis[i] <- sum(tab$D.fis > (0.75)) / 100 
    ptab$pvalue.fus[i] <- sum(tab$D.fus > (0.75)) / 100 
    ptab$pvalue.comb[i] <- sum(tab$D.comb > (0.75)) / 100
  }
  if(ptab$condition[i] %in% c(2,4,6) & ptab$rr[i] == 5){
    ptab$pvalue.fis[i] <- sum(tab$D.fis > (0.75 * 4)) / 100 
    ptab$pvalue.fus[i] <- sum(tab$D.fus > (0.75 * 4)) / 100 
    ptab$pvalue.comb[i] <- sum(tab$D.comb > (0.75 * 4)) / 100
  }
  if(ptab$condition[i] %in% c(2,4,6) & ptab$rr[i] == 10){
    ptab$pvalue.fis[i] <- sum(tab$D.fis > (0.75 * 9)) / 100 
    ptab$pvalue.fus[i] <- sum(tab$D.fus > (0.75 * 9)) / 100
    ptab$pvalue.comb[i] <- sum(tab$D.comb > (0.75 * 9)) / 100
  }
  if(ptab$condition[i] %in% c(1,3,5)){
    ptab$pvalue.fis[i] <- sum(tab$D.fis > (0)) / 100 
    ptab$pvalue.fus[i] <- sum(tab$D.fus > (0)) / 100
    ptab$pvalue.comb[i] <- sum(tab$D.comb > (0)) / 100
  }
  ptab$passZeroFis[i] <- sum(tab$D.fis.zero) / 100
  ptab$passZeroFus[i] <- sum(tab$D.fus.zero) / 100
  ptab$passZeroComb[i] <- sum(tab$D.comb.zero) / 100
}
# re order ptab by tip number
ptab <- ptab[order(ptab$ntip),]
# rename conditions 5 and 6 tip values accordingly
ptab$ntip[ptab$condition %in% c(5,6)] <- rep(c(2,4,6,8,10), each = 6)
#### plot ####
par(mfrow = c(2,3))
#set colours and point types
col <- viridis::viridis(4, alpha = 0.5)

#### plot fission ####
## empirical p Values ##
### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 2") 
for(i in 1:4){
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$pvalue.fis[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
  
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

# add legend
legend("topleft",inset=.02,
       legend=c("Condition 1",
                "Condition 2",
                "Condition 3",
                "Condition 4"),
       col=col, 
       lty=1,
       cex=1,
       lwd = 2)
### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,.6),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 5") 
for(i in 1:4){
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$pvalue.fus[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 10") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$pvalue.fus[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

## Proportion of positive resultss ##
### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$passZeroFis[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
  
}
### rates ratio 5 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips",
     main = expression(paste(Delta, "R"[Fission]))) 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$passZeroFis[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
  
}
### rates ratio 10 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$passZeroFis[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}

#### condition 5 and 6 fission #### 
par(mfrow = c(2,3))
## emperical p values ##
## rates ratio 2 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 2")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$pvalue.fis[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$pvalue.fis[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

# add legend
legend("topleft",
       inset=.02,
       legend=c("Condition 5",
                "Condition 6"),
       col=c("#f1a340","#998ec3"),
       lty=1,
       cex=1,
       lwd = 2)

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 5")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$pvalue.fis[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$pvalue.fis[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 10")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$pvalue.fis[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$pvalue.fis[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

## Proportion of positive resultss ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$passZeroFis[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$passZeroFis[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates",
     main = expression(paste(Delta, "R"[Fission])))

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$passZeroFis[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$passZeroFis[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$passZeroFis[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$passZeroFis[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
#### plot fusion ####
## empirical p Values ##
### rates ratio 2 ###
par(mfrow = c(2,3))
### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 2") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$pvalue.fis[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add legend
legend("topleft",inset=.02,
       legend=c("Condition 1",
                "Condition 2",
                "Condition 3",
                "Condition 4"),
       col=col, 
       lty=1,
       cex=1,
       lwd = 2)
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
### rates ratio 5 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 5") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$pvalue.fis[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
### rates ratio 10 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 10") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$pvalue.fis[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## Proportion of positive resultss ##
### rates ratio 2 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 

for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$passZeroFus[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
### rates ratio 5 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips",
     main = expression(paste(Delta, "R"[Fusion]))) 

for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$passZeroFus[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
### rates ratio 10 ###
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 

for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$passZeroFus[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}

#### condition 5 and 6 fusion #### 
par(mfrow = c(2,3))
## emperical p values ##
## rates ratio 2 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 2")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$pvalue.fus[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$pvalue.fus[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
# add legend
legend("topleft",
       inset=.02,
       legend=c("Condition 5",
                "Condition 6"),
       col=c("#f1a340","#998ec3"),
       lty=1,
       cex=1,
       lwd = 2)

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 5")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$pvalue.fus[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$pvalue.fus[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 10")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$pvalue.fus[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$pvalue.fus[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## Proportion of positive resultss ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates",)

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$passZeroFus[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$passZeroFus[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates",
     main = expression(paste(Delta, "R"[Fusion])))

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$passZeroFus[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$passZeroFus[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$passZeroFus[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$passZeroFus[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

#### plot combined  conditions 1 to 4 ####
# plot
par(mfrow = c(2,3))
## empirical p Values
## RR 2
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 2") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$pvalue.comb[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
  
}
# add legend
legend("topleft",inset=.02,
       legend=c("Condition 1",
                "Condition 2",
                "Condition 3",
                "Condition 4"),
       col=col, lty=1, cex = 1,
       lwd = 2)
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## RR 5
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 5") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$pvalue.comb[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## RR 10
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Number of tips",
     main = "Rates ratio = 10") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$pvalue.comb[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## Proportion of positive resultss
## RR 2
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 2],
        y = ptab$passZeroComb[ptab$condition == i & ptab$rr == 2],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
  
}
## RR 5
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips",
     main = expression(paste(Delta, "R"[combined]))) 

for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 5],
        y = ptab$passZeroComb[ptab$condition == i & ptab$rr == 5],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}
## RR 10
plot(x = NULL,
     y = NULL,
     xlim = c(50,250),
     ylim = c(0,1),
     ylab = "Proportion of positive results",
     xlab = "Number of tips") 
for(i in 1:4){
  # plot condition 2
  lines(x = ptab$ntip[ptab$condition == i & ptab$rr == 10],
        y = ptab$passZeroComb[ptab$condition == i & ptab$rr == 10],
        pch = 16,
        type = "o", lwd = 2,
        col = col[i])
}


#### condition 5 and 6 combined #### 
par(mfrow = c(2,3))
## emperical p values ##
## rates ratio 2 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 2")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$pvalue.comb[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$pvalue.comb[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

# add legend
legend("topleft",
       inset=.02,
       legend=c("Condition 5",
                "Condition 6"),
       col=c("#f1a340","#998ec3"),
       lty=1,
       cex=1,
       lwd = 2)
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 5")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$pvalue.comb[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$pvalue.comb[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)
## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Emperical p-value",
     xlab = "Percentage of tips \nhaving differential rates",
     main = "Rates ratio = 10")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$pvalue.comb[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$pvalue.comb[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")
# add 0.05 p value cut off
abline(h = 0.05, col = "red", lty = 2)

## Proportion of positive resultss ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 2],
      y = ptab$passZeroComb[ptab$condition == 5 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 2],
      y = ptab$passZeroComb[ptab$condition == 6 & ptab$rr == 2],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 5 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates",
     main = expression(paste(Delta, "R"[combined])))

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 5],
      y = ptab$passZeroComb[ptab$condition == 5 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 5],
      y = ptab$passZeroComb[ptab$condition == 6 & ptab$rr == 5],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

## rates ratio 10 ##
plot(x = NULL,
     y = NULL,
     xlim = c(2,10),
     ylim = c(0,1),
     ylab = "Proportion of positive resultss",
     xlab = "Percentage of tips \nhaving differential rates")

# plot condition 5
lines(x = ptab$ntip[ptab$condition == 5 & ptab$rr == 10],
      y = ptab$passZeroComb[ptab$condition == 5 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#f1a340")

# plot condition 6
lines(x = ptab$ntip[ptab$condition == 6 & ptab$rr == 10],
      y = ptab$passZeroComb[ptab$condition == 6 & ptab$rr == 10],
      pch = 16,
      type = "o", lwd = 2,
      col = "#998ec3")

