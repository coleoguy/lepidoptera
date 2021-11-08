# Terrence Sylvester
# 28 August 2021
# pradakshanas@gmail.com

#load libraries
library(ape)
library(phytools)
library(plotrix)
library(ggplot2)

# helper functions
source("../rscripts/helper.functions.R")

#load data
load("../results/7a.Lepidoptera.branch.specific.rates.RData")

# make three data tables
tipRates <- tipChroms <- tipFeeding <- as.data.frame(matrix(data = NA,
                                                            ncol = 100,
                                                            nrow = Ntip(tree)))
# get colnames
colnames(tipRates) <- colnames(tipChroms) <- colnames(tipFeeding) <- paste("test",1:100)
# get row names
row.names(tipRates) <- row.names(tipChroms) <- row.names(tipFeeding) <- tree$tip.label
# fill each table

for(j in 1:100){
  tipRates[,j] <- final.results[[j]][row.names(tipRates),5]
  tipChroms[,j] <- final.results[[j]][row.names(tipChroms),2]
  tipFeeding[,j] <- final.results[[j]][row.names(tipFeeding),4]
}

# get the final result from all 100 runs
tipRateStat <- as.data.frame(matrix(data = NA,
                                    ncol = 5,
                                    nrow = Ntip(tree)))
colnames(tipRateStat) <- c("species", "genProb","chromNum","sdChrom","rate")
# fill
tipRateStat$species <- row.names(tipRates)
tipRateStat$genProb <- rowMeans(tipFeeding)
tipRateStat$rate <- rowMeans(tipRates)
tipRateStat$chromNum <- rowMeans(tipChroms)

# find species and rates with mixed binary trait
tipRateStat[tipRateStat$genProb > 0 & tipRateStat$genProb < 1,]

# assign a value to species that have a mixed binary trait
tipRateStat$genProb[tipRateStat$genProb > 0 & tipRateStat$genProb < 1] <- 2

for(i in 1:441){
  tipRateStat$sdChrom[i] <- sd(as.numeric(tipChroms[i,]))
}

# plot tree with branches
scale <- .9
plotTree.wBars(tree = tree,
               x = setNames(tipRateStat$rate + 2, tipRateStat$species),
               type = "fan",
               col = viridis::viridis(3, end = 0.8, alpha = .7)[tipRateStat$genProb+1],
               border = NA,
               lwd = 1,
               scale = scale,
               width = 1)
radi <- getRadius(scale = scale,
                  width = 1,
                  tree = tree,
                  tip.labels = FALSE,
                  trait.values = tipRateStat$rate + 2,
                  classes = 3)
# get radius for zero line
radi_0 <-  getRadius(scale = scale,
                     width = 1,
                     tree = tree,
                     tip.labels = FALSE,
                     trait.values = 2,
                     classes = 1)

#plot circles
draw.circle(x = 0,
            y = 0,
            radius = radi,
            nv=100,
            border="gray",
            col=NA,
            lty=2,
            density=NULL,
            angle=45,
            lwd=1)

# draw circle
draw.circle(x = 0,
            y = 0,
            radius = radi_0,
            nv=100,
            border="red",
            col=NA,
            lty=2,
            density=NULL,
            angle=45,
            lwd=1)

# add corresponding rate class for each rate
text(x = radi,
     y = rep(0, length(radi)),
     labels = round(as.numeric(names(radi))-2,0),
     srt = -90, 
     cex = 0.7,
     pos = 4)

text(x = radi_0,
     y = 0,
     labels = 0,
     srt = -90, 
     cex = 0.7,
     pos = 4)


# add legend
xPoints <- max(radi)

off <- 0.1

text(x = xPoints - (xPoints* off),
     y = xPoints - (xPoints* off),
     labels = "Binary trait",
     pos = 4,
     cex = .9)

points(x = rep(xPoints - (xPoints* off),3),
       y = c((xPoints - (xPoints* off)) - xPoints * 0.05,
             (xPoints - (xPoints* off)) - xPoints * 0.1,
             (xPoints - (xPoints* off)) - xPoints * 0.15),
       pch = 16,
       col = c(viridis::viridis(3,end = 0.8,direction = -1)[[2]],
               viridis::viridis(3,end = 0.8,direction = -1)[[1]],
               viridis::viridis(3,end = 0.8,direction = -1)[[3]]))
text(x = rep(xPoints - (xPoints* off),3),
     y = c((xPoints - (xPoints* off)) - xPoints * 0.05,
           (xPoints - (xPoints* off)) - xPoints * 0.1,
           (xPoints - (xPoints* off)) - xPoints * 0.15),
     labels = c("Generalists","Mixed", "Specialists"),
     pos = 4,
     cex = 0.8)

# non parametric t test
# all rates
wilcox.test(tipRateStat$rate ~ tipRateStat$genProb)
# highest rate removed
wilcox.test(tipRateStat$rate[-c(133)] ~ tipRateStat$genProb[-c(133)])
# top 3 rates removed
wilcox.test(tipRateStat$rate[-c(133,404,299)] ~ tipRateStat$genProb[-c(133,404,299)])
# top 5% removed
wilcox.test(tipRateStat$rate[(tipRateStat$rate <= quantile(tipRateStat$rate, 0.95))] ~ tipRateStat$genProb[(tipRateStat$rate <= quantile(tipRateStat$rate, 0.95))])
# only in top 5%
wilcox.test(tipRateStat$rate[(tipRateStat$rate > quantile(tipRateStat$rate, 0.95))] ~ tipRateStat$genProb[(tipRateStat$rate > quantile(tipRateStat$rate, 0.95))])

# log transformed dodged bar plot
ggplot(tipRateStat, aes(x = (rate)^2, 
                        fill = as.factor(genProb))) +
  geom_histogram(position = "dodge2",
                 binwidth = 25) +
  scale_fill_discrete("Larval feeding type",
                      labels = c("Specialists", "Generalists"), 
                      type = c("#f1a340","#998ec3")) +
  xlab(expression(paste("Species rate ", "(MY"^-1,")"))) + 
  ylab("Count") + 
  theme_bw()


# log transformed woth zero rates dodged bar plot
ggplot(tipRateStat, aes(x = log(rate + 1), 
                        fill = as.factor(genProb))) +
  geom_histogram(position = "dodge2",
                 binwidth = 0.25) +
  scale_fill_discrete("Larval feeding type",
                      labels = c("Specialists", "Generalists"), 
                      type = c("#f1a340","#998ec3")) +
  xlab(expression(paste("Species rate ", "(MY "^-1,")"))) + 
  ylab("Count") + 
  scale_x_continuous(breaks = log(c(0.01,0.1,1,10,100,1000)), 
                     labels = c(0.01,0.1,1,10,100,1000)) +
  theme_bw()

# square root transformed dodged bar plot
ggplot(tipRateStat, aes(x = sqrt(rate), 
                        fill = as.factor(genProb))) +
  geom_histogram(position = "dodge2",
                 binwidth = 0.25) +
  scale_fill_discrete("Larval feeding type",
                      labels = c("Specialists", "Generalists"), 
                      type = c("#f1a340","#998ec3")) +
  xlab(expression(paste("Species rate ", "(MY "^ -1,")"))) + 
  ylab("Count") + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12), 
                     labels = c(0,1,4,9,16,25,36,49,64,81,100,121,144)) +
  theme_bw()
