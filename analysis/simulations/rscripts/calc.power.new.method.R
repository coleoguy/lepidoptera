# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# load libraries
library(coda)
library(stringr)
library(doSNOW)

# get helper functions
source("../rscripts/functions.R")

# set pars for multitasking
# define number of clusters for parallel computing
NumberOfClusters <- 4
cl <- makeCluster(NumberOfClusters, outfile = "")
registerDoSNOW(cl)
# get empirical p values
nrep <- length(dir("../results/empP-MCMC/"))
# get the working tree number
wtree <- dir("../results/empP-MCMC/rep-1/")
wtree <- as.numeric(gsub("tree-", "", wtree))
# get number of files in a single repetition
nfiles <- length(dir(paste("../results/empP-MCMC/rep-", sample(nrep,1), "/tree-",
                           wtree, sep = "")))
x <- foreach(i = 1:nrep, .verbose = T,.packages = c("stringr"))%dopar%{
  files <- dir(paste("../results/empP-MCMC/rep-", i, "/tree-", wtree, sep = "")) 
  # make the data table
  dat <- as.data.frame(matrix(data = NA, nrow = nfiles, ncol = 8))
  colnames(dat) <- c("run","rep","fission","fusion","aneuploidy", "tips","rr","cond")
  for(ii in 1:nfiles){
    # load results
    load(paste("../results/empP-MCMC/rep-", i, "/tree-", wtree, "/", files[ii],
               sep = ""))
    dat$run[ii] <- files[ii]
    dat$rep[ii] <- i
    dat$fission[ii] <- empP$empiricalP$EmpPvalueFission
    dat$fusion[ii] <- empP$empiricalP$EmpPvalueFusion
    dat$aneuploidy[ii] <- empP$empiricalP$EmpPvalueAneuploidy
    dat$tips[ii] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[ii], "([A-z])", replacement = ""), split = ".",fixed = T))[3])
    dat$rr[ii] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[ii], "([A-z])", replacement = ""), split = ".", fixed = T))[4])
    dat$cond[ii] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[ii], "([A-z])", replacement = ""),split = ".", fixed = T))[2])
    # remove unwanted results
    rm(empP)
  }
  dat
}
# stop cluster
stopCluster(cl)
# combine results from multi threading into a single data frame
dat <- x[[1]]
for(i in 2:length(x)){
  dat <- rbind(dat,x[[i]])
}
# sort data table by tip number
dat <- dat[order(dat$tips),]
# calculate the number of times we get a significant result
empPower <- as.data.frame(matrix(data = NA, nrow = nfiles, ncol = 7))
colnames(empPower) <-  c("run", "ntip", "cond", "rr", "fission","fusion", "aneuploidy")
files <- dir(paste("../results/empP-MCMC/rep-", sample(nrep,1), "/tree-", wtree, sep = "")) 
for(i in 1:nrow(empPower)){
  empPower$run[i] <- files[i]
  empPower$ntip[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i], "([A-z])", replacement = ""), split = ".", fixed = T))[3])
  empPower$cond[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i],"([A-z])", replacement = ""), split = ".", fixed = T))[2])
  empPower$rr[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i], "([A-z])", replacement = ""), split = ".", fixed = T))[4])
  empPower$fusion[i] <- sum(dat$fission[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
  empPower$fission[i] <- sum(dat$fusion[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
  empPower$aneuploidy[i] <- sum(dat$aneuploidy[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
}
# save data
write.csv(dat,"../results/power.new.method.full.tab.23-rep.csv", row.names = F)
write.csv(empPower, "../results/power.new.method.23-rep.csv",row.names = F)
