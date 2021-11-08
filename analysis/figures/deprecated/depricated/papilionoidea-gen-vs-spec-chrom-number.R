# Terrence Sylvester
# pradakshanas@gmail.com
# May 18, 2020
# distribution of chromosome numbers

# load libraries
library(ape)

# load helper functions
source("../rscripts/helper.functions.R")

# read in chromosome data
chroms <- read.csv("../data/chroms/lep.chroms.csv", 
                   as.is = T)

# read in tree data
trees <- read.tree("../data/trees/papilionoidea-10.trees")

# read in hosts data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt",
                  as.is = T)

# sampe chromosome counts
chroms <- chromSampler(dat = chroms)

# get the trait overlap
dat <- TraitOverlap(dat = chroms,
                    trees = trees,
                    hosts = hosts)

# make chromosome number haploid
dat$haploid <- dat$haploid / 2


# remove all but keep dat
rm(list = ls()[-5])

# in this data set specialists are 1 and generalists are 0

# plot the distribution
# make the plot region
plot(x = NULL,
     y = NULL,
     xlim = c(0,1),
     ylim = c(0, max(dat$haploid)),
     xlab = "Host type",
     ylab = "Haploid chromosome count",xaxs = "i",
     xaxt = "n")

# plot tick marks
axis(side = 1,
     at= c(0.2,0.8), 
     labels= c("Generalist species", "Specialist species"))

# plot generalists
points(x = jitter(rep(.2, sum(dat$hosts == 0)), 12),
       y = dat$haploid[dat$hosts == 0],
       pch = 16,
       col = rgb(0,0,0,.5))

# plot specialists
points(x = jitter(rep(.8, sum(dat$hosts == 1)), 3),
       y = dat$haploid[dat$hosts == 1],
       pch = 16,
       col = rgb(0,0,0,.5))


