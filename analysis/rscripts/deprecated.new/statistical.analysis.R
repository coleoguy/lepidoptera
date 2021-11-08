# Terrence Sylvester
# pradakshanas@gmail.com
# 21 May 2020
# ANOVA

# load libraries
library(phytools)

# load helper functions
source("../rscripts/helper.functions.R")

# read in data
# read chromosome data
chroms <- read.csv("../data/chroms/lep.chroms.csv", as.is = T)

# read trees
trees <- read.tree("../data/trees/papilionoidea-10.trees")

# read host data
hosts <- read.csv("../data/hosts/papilionoidea-bisse-data.txt", as.is = T)

# sample chromosomes
dat <- chromSampler(dat = chroms)

# trait overlap
dat <- TraitOverlap(dat = dat, trees = trees, hosts = hosts)

# make sure haploid corrected 
dat$haploid <- dat$haploid / 2

# chainge tip names
tree <- changeTipnames(tree = trees[[sample(1:10, 1)]], dat = dat)

## ANOVA 
aov.out <- aov(dat$haploid ~ as.factor(dat$hosts), data = dat)
summary(aov.out)

# post hoc test to see which has higher variance
TukeyHSD(aov.out)

## phylo anova
groups <- dat$hosts
response <- dat$haploid
names(groups) <- names(response) <- dat$species

phy.aov.out <-  phylANOVA(tree = tree,
                          x = groups,
                          y = response,
                          p.adj = "bonferroni")

phy.aov.out

## t-test to see the differnce in the means
t.test.out <- t.test(dat$haploid[dat$hosts == 0], dat$haploid[dat$hosts==1])
t.test.out