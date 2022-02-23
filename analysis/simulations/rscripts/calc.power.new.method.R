# Terrence Sylvester
# pradakshanas@gmail.com
# 3rd February 2022

# load libraries
library(coda)
library(stringr)
# get helper functions
source("../rscripts/functions.R")

# get empirical p values
nrep <- length(dir("../results/empP-MCMC/"))
# get the working tree number
wtree <- dir("../results/empP-MCMC/rep-1/")
wtree <- as.numeric(gsub("tree-", "", wtree))
# get number of files in a single repetition
nfiles <- length(dir(paste("../results/empP-MCMC/rep-",
                           sample(nrep,1),
                           "/tree-",
                           wtree, 
                           sep = "")))
# make the data table
dat <- as.data.frame(matrix(data = NA,
                            nrow = nrep * nfiles,
                            ncol = 8))
colnames(dat) <- c("run","rep","fission","fusion","aneuploidy", "tips","rr","cond")
# fill in data table
counter <- 1
for(i in 1:nrep){
  files <- dir(paste("../results/empP-MCMC/rep-",
                     i,
                     "/tree-",
                     wtree, 
                     sep = "")) 
  for(ii in 1:nfiles){
    cat("working on file",counter, "out of", nrow(dat), "\n")
    # load results
    load(paste("../results/empP-MCMC/rep-",
               i,
               "/tree-",
               wtree,
               "/",
               files[ii],
               sep = ""))
    dat$run[counter] <- files[ii]
    dat$rep[counter] <- i
    dat$fission[counter] <- empP$empiricalP$EmpPvalueFission
    dat$fusion[counter] <- empP$empiricalP$EmpPvalueFusion
    dat$aneuploidy[counter] <- empP$empiricalP$EmpPvalueAneuploidy
    dat$tips[counter] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[counter],
                                                                    "([A-z])",
                                                                    replacement = ""),
                                                    split = ".",
                                                    fixed = T))[3])
    dat$rr[counter] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[counter],
                                                                  "([A-z])",
                                                                  replacement = ""),
                                                  split = ".",
                                                  fixed = T))[4])
    dat$cond[counter] <- as.numeric(unlist(strsplit(str_replace_all(dat$run[counter],
                                                                    "([A-z])",
                                                                    replacement = ""),
                                                    split = ".",
                                                    
                                                    fixed = T))[2])
    counter <- counter + 1
    print(counter)
    # remove unwanted results
    rm(empP)
  }
}
# sort data table by tip number
dat <- dat[order(dat$tips),]

# calculate the number of times we get a significant result
empPower <- as.data.frame(matrix(data = NA,
                                 nrow = nfiles,
                                 ncol = 7))
colnames(empPower) <-  c("run", "ntip", "cond", "rr", "fission","fusion", "aneuploidy")

for(i in 1:nrow(empPower)){
  empPower$run[i] <- files[i]
  
  empPower$ntip[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i],
                                                                 "([A-z])",
                                                                 replacement = ""),
                                                 split = ".",
                                                 fixed = T))[3])
  empPower$cond[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i],
                                                                 "([A-z])",
                                                                 replacement = ""),
                                                 split = ".",
                                                 fixed = T))[2])
  empPower$rr[i] <- as.numeric(unlist(strsplit(str_replace_all(empPower$run[i],
                                                               "([A-z])",
                                                               replacement = ""),
                                               split = ".",
                                               fixed = T))[4])
  
  empPower$fusion[i] <- sum(dat$fission[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
  empPower$fission[i] <- sum(dat$fusion[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
  empPower$aneuploidy[i] <- sum(dat$aneuploidy[dat$tips == empPower$ntip[i] &  dat$cond == empPower$cond[i] & dat$rr == empPower$rr[i]] <= 0.05) / nrep
}
# save data
write.csv(dat,"../results/power.new.method.full.tab.csv", row.names = F)
write.csv(empPower, "../results/power.new.method.csv",row.names = F)


