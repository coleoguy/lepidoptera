# Terrence Sylvester
# pradakshanas@gmail.com
# May 07, 2020
# mcmc analysis - burnin

# load libraries
library(viridis)

# load data
load("../results/papilionoidea-mcmc.run-Mkn.RData")

# remove unwanted data
rm(list = ls()[-14])

# for the time being we will not concider the results from tree two as this 
# was a failed mcmc run

# look at the mcmc
plot(results[[1]]$p, 
     type = "l",
     ylim = c(-820,-750))

for(i in 2:10){
  lines(results[[i]]$p,
        col = viridis(n = 10)[i])
}

# clear plot window
dev.off()

# looking at this we can see that the mcmc runs rached convergence by 200
# generations. Therefore we will discard initial 25% as burnin and sample
# from the remaining 75%

# define burnin
burn.in <- .25
pre.burn.in <- seq(from = 1,
                   to = nrow(results[[1]]) * burn.in,
                   by = 1)

# remove burnin and sample 10000 data points from post burnin
post.burn.in <- as.data.frame(matrix(data = NA, 
                                     nrow = nrow(results[[1]]) * (1-burn.in) * 10,
                                     ncol = ncol(results[[1]])))

colnames(post.burn.in) <- colnames(results[[1]])


for(i in 1:10){
  x <- results[[i]][-pre.burn.in,]
  rng <- which(is.na(post.burn.in[,1]))[1]:(which(is.na(post.burn.in[,1]))[1] + nrow(x) - 1)
  post.burn.in[(rng),] <- x 
}

# remove empty rows from post burn in
# post.burn.in <- post.burn.in[-which(is.na(post.burn.in$i)),]

# clear all but post burn in
rm(list = ls()[-3])

# save results as a csv file
write.csv(x = post.burn.in,
          file = "../results/papilionoidea-mcmc.run-Mkn-post-burn-in.csv",
          row.names = F)
