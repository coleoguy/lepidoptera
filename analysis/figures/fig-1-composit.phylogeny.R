# Terrence Sylvester
# 28 August 2021
# pradakshanas@gmail.com

# load required libraries
library(ape)
library(phytools)

# load helper functions
source("../rscripts/helper.functions.R")

# read in data
trait <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
chrom <- read.delim("../data/chroms/chroms.txt", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# fill in binomial column so that it does not includes sub species names
for(i in 1:nrow(chrom)){
        chrom$binomial[i] <- paste(chrom$Genus[i], chrom$species[i], sep = "_")
}

# filter out unwanted columns
chrom <- chrom[c("Family","Genus","binomial","female2N", "male2N", "notesMale2N")]

# process data
# sample a single chromosome number when there is a range
proc.chrom <- chromSampler(dat = chrom)
# remove species that lack male chromosome number data
proc.chrom <-  proc.chrom[!(is.na(proc.chrom$male2N)),]
# remove unwanted columns from the data table
proc.chrom <- proc.chrom[c("Family","Genus","binomial","male2N")]
# remove species that do not have information on chromosome number
proc.chrom <- proc.chrom[!is.na(proc.chrom$male2N),]
# get the haploid chromosome number
proc.chrom$haploid <- NA
proc.chrom$haploid <- proc.chrom$male2N / 2

# if there are species with floating point values for haploid number. 
# randomly sample a single value (ceiling or floor) for these species
for(i in 1:nrow(proc.chrom)){
        if(proc.chrom$haploid[i] %% 1 != 0){
                proc.chrom$haploid[i] <- sample(c(ceiling(proc.chrom$haploid[i]), floor(proc.chrom$haploid[i])),1)
        }
}
# check if there are species with multiple records for chromosome number
# if present, randomly sample a single value
chrom.dedup <- as.data.frame(matrix(data = NA,
                                    nrow = length(unique(proc.chrom$binomial)),
                                    ncol = ncol(proc.chrom)))
colnames(chrom.dedup) <- colnames(proc.chrom)

# fill out species names
chrom.dedup$binomial <- unique(proc.chrom$binomial)
# use a for loop to fill out rest of the information
# for species that have multiple records we first store them in a temporary object
# and then sample a single value
for (i in 1:nrow(chrom.dedup)) {
        if(sum(proc.chrom$binomial == chrom.dedup$binomial[i]) == 1){
                chrom.dedup$Family[i] <- proc.chrom$Family[proc.chrom$binomial == chrom.dedup$binomial[i]]
                chrom.dedup$Genus[i] <- proc.chrom$Genus[proc.chrom$binomial == chrom.dedup$binomial[i]]
                chrom.dedup$haploid[i] <- proc.chrom$haploid[proc.chrom$binomial == chrom.dedup$binomial[i]]
                
        }
        if(sum(proc.chrom$binomial == chrom.dedup$binomial[i]) > 1){
                hit <- proc.chrom$haploid[proc.chrom$binomial == chrom.dedup$binomial[i]] 
                chrom.dedup$haploid[i] <- sample(hit, 1)
        }
        chrom.dedup$Family[i] <- sample(proc.chrom$Family[proc.chrom$binomial == chrom.dedup$binomial[i]],1)
        chrom.dedup$Genus[i] <- sample(proc.chrom$Genus[proc.chrom$binomial == chrom.dedup$binomial[i]],1)
}
# remove unwanted columns
chrom.dedup <- chrom.dedup[,-4]
# get the development state of each species
chrom.dedup$feeding <- NA
for(i in 1:nrow(chrom.dedup)){
        hit.feeding <- trait$type[trait$binomial == chrom.dedup$binomial[i]]
        if(length(hit.feeding) != 0){
                chrom.dedup$feeding[i] <- trait$type[trait$binomial == chrom.dedup$binomial[i]]
        }
}
# get the genera level representations for development
# get genera which have no development data
gen.feeding.no <- chrom.dedup$Genus[is.na(chrom.dedup$feeding)]
# get genera which have development data
gen.feeding.yes <- chrom.dedup$Genus[!is.na(chrom.dedup$feeding)]
# get genera that are absent the subset of genera that have development data
gen.feeding.no.match <- c()
for(i in 1:length(gen.feeding.no)){
        if(gen.feeding.no[i] %in% gen.feeding.yes == F){
                gen.feeding.no.match <- c(gen.feeding.no.match, gen.feeding.no[i])
        }
}
# get the unique list of names of these genera
gen.feeding.no.match.unique <- unique(gen.feeding.no.match)
# fill these genera with development data
for(i in 1:length(gen.feeding.no.match.unique)){
        if(gen.feeding.no.match.unique[i] %in% chrom.dedup$Genus){
                hit.feeding.genera <- trait$type[trait$genus == gen.feeding.no.match.unique[i]]
                if(length(hit.feeding.genera) != 0){
                        if(length(unique(hit.feeding.genera)) == 1){
                                chrom.dedup$feeding[chrom.dedup$Genus == gen.feeding.no.match.unique[i]] <- unique(hit.feeding.genera)
                        }else{
                                chrom.dedup$feeding[chrom.dedup$Genus == gen.feeding.no.match.unique[i]] <- sample(hit.feeding.genera,1)
                        }
                }
        }
}
# remove species that do not have Development information
chrom.dedup <- chrom.dedup[!is.na(chrom.dedup$feeding),]
# get the overlap between the species and chromosome number data sets
sp.overlap <-  sp.matches(dat = chrom.dedup, trees = trees, hyper = T)
# processed dat
processed.dat <- sp.overlap$chroms
# rename tree tip names to reflect genera level matches
tree <- trees[[sample(1:length(trees),1)]]
tree <- keep.tip(tree, sp.overlap$name_corrections$name_on_tree)
for(i in 1:Ntip(tree)){
        tree$tip.label[i] <- sp.overlap$name_corrections$species_name[sp.overlap$name_corrections$name_on_tree == tree$tip.label[i]]
}
# sort the chromosome number dataset so that it matches with the tip order
chrom.sorted <- as.data.frame(matrix(data = NA,
                                     nrow = nrow(processed.dat),
                                     ncol = ncol(processed.dat)))
colnames(chrom.sorted) <- colnames(processed.dat)
for(i in 1:nrow(chrom.sorted)){
        chrom.sorted[i,] <- processed.dat[which(processed.dat$binomial == tree$tip.label[i]),]
}
# make a new data table which has species names and chromosome number
MCMC.dat <- as.data.frame(matrix(data = NA,
                                 nrow = nrow(chrom.sorted),
                                 ncol = 3))
colnames(MCMC.dat) <- c("Species", "Chroms","feeding")
MCMC.dat$Species <- chrom.sorted$binomial
MCMC.dat$Chroms <- chrom.sorted$haploid
MCMC.dat$feeding <- chrom.sorted$feeding
# in the dataset of larval feeding 1 means specialists and 0 means generalists
# change that so 1 means generalist and 0 means specialist
MCMC.dat$gen.prob <- NA
for(i in 1:nrow(MCMC.dat)){
        if(MCMC.dat$feeding[i] == 1){
                MCMC.dat$gen.prob[i] <- 0
        }
        if(MCMC.dat$feeding[i] == 0){
                MCMC.dat$gen.prob[i] <- 1
        }
}
# # plot the tree
# plotTree.wBars(tree = tree,
#                x = setNames(MCMC.dat$Chroms, MCMC.dat$Species),
#                type = "fan",
#                scale = .2,
#                method = "plotTree",
#                width = 1,
#                offset = 1,
#                border = NA,
#                lwd = 1)

# get the contmap of chromosome number
map <- contMap(tree = tree,
               x = setNames(log(MCMC.dat$Chroms), MCMC.dat$Species),
               res = 1000,
               ftype = "off",
               plot = F)
#change the colour scale of the contMap
map$cols <- viridis::viridis(n = length(map$cols))
names(map$cols) <- 0:(length(map$cols) -1)

# plot
plot.contMap(map, 
             type = "fan", 
             ftype = "off",
             outline = F,
             lwd = 1.5,
             legend = F,
             res = 1000,
             mar = rep(.6, 4))
add.color.bar(50,
              map$cols,
              title="",
              lims=range(MCMC.dat$Chroms),
              digits=3,
              prompt=FALSE,
              x=-150, # this is the possition where legend starts. We chose this
              # by using the locator function in R
              y=-80,
              lwd=4,
              fsize=1,
              subtitle="Chromosome number")

# get the overlap between the species and chromosome number data sets
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
n<-Ntip(tree)
col.bin<-setNames(c("red","blue"),levels(as.factor(MCMC.dat$gen.prob)))

# mark tips
offset <- 0.03
par(lend=3)
for(i in 1:Ntip(tree)){
        if(max(branching.times(tree)) <= 1){
                cc<-if(obj$xx[i]>0) 5 else -5
                th<-atan(obj$yy[i]/obj$xx[i])
                segments(obj$xx[i] + (obj$xx[i] * offset),
                         obj$yy[i] + (obj$yy[i] * offset),
                         obj$xx[i]+(cc/tree.depth)*cos(th) + (obj$xx[i] * offset),
                         obj$yy[i]+(cc/tree.depth)*sin(th) + (obj$yy[i] * offset),
                         lwd=2,
                         col=col.bin[MCMC.dat$gen.prob[MCMC.dat$Species == tree$tip.label[i]] + 1])
                
        }
        else{
                cc<-if(obj$xx[i]>0) 5 else -5
                th<-atan(obj$yy[i]/obj$xx[i])
                segments(obj$xx[i] + (obj$xx[i] * offset),
                         obj$yy[i] + (obj$yy[i] * offset),
                         obj$xx[i]+cc*cos(th) + (obj$xx[i] * offset),
                         obj$yy[i]+cc*sin(th) + (obj$yy[i] * offset),
                         lwd=2,
                         col=col.bin[MCMC.dat$gen.prob[MCMC.dat$Species == tree$tip.label[i]] + 1])
        }
}

# legend
# Binary state
text(x = 95, y = 90, labels = "Binary State",
     pos = 4)
points(x = rep(95, 2),
       y = c(84,78),
       pch = 16,
       # bg = c("#4daf4a","#984ea3"),
       col = c("blue","red"),
       cex = 1)
text(x = rep(95, 2),
     y = c(84,78),
     labels = c("Generalists",
                "Specialits"),
     pos = 4,
     cex = .8)

