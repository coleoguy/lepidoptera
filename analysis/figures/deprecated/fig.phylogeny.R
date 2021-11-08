# Terrence Sylvester
# 14 Oct 2020
# pradakshanas@gmail.com

# load libraries
library(phytools)
library(viridis)
library(chromePlus)
library(diversitree)

# get helper functions
source("../rscripts/helper.functions - Copy.R")

# read in data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")

# combine genus and species names to get the full name in hosts dataset and in
# chroms dataset
host$SpeciesNames <- paste(host$genus, host$species, sep = "_")

# fill in species for data
dat$SpecisName <- paste(dat$other.names.genera,
                        dat$other.names.species,
                        sep = "_")
dat$SpecisName[dat$SpecisName == "_"] <- paste(dat$Genus[dat$SpecisName == "_"],
                                               dat$species[dat$SpecisName == "_"],
                                               sep = "_")


# sample a single chromosome number when there is a range
new.dat <- chromSampler(dat = dat)
# remove species that lack male chromosome number data
new.dat <-  new.dat[!(is.na(new.dat$male2N)),]
# get the overlap between the traits and the phylogeny
overlap <- sp.matches(dat = new.dat,
                      trees = trees)
# proccessed dat
processed.dat <- overlap$chroms
# rename tree and hosts species names
tree <- trees[[1]]
tree <- keep.tip(tree, overlap$name_corrections$name_on_tree)
# rename tree tip and host species name to reflect genera level matches
for(i in 1:Ntip(tree)){
  tree$tip.label[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == tree$tip.label[i]]
}
host <- host[host$SpeciesNames %in% overlap$name_corrections$name_on_tree,]
for(i in 1:nrow(host)){
  host$SpeciesNames[i] <- overlap$name_corrections$species_name[overlap$name_corrections$name_on_tree == host$SpeciesNames[i]]
}
# make a new data table which has species names chromosome number and host type
MCMC.dat <- as.data.frame(matrix(data = NA,
                                 nrow = nrow(processed.dat),
                                 ncol = 3))
colnames(MCMC.dat) <- c("Species", "Chroms", "HostType")
MCMC.dat$Species <- processed.dat$SpecisName
MCMC.dat$Chroms <- processed.dat$male2N / 2
for(i in 1:nrow(MCMC.dat)){
  MCMC.dat$HostType[i] <- host$type[host$SpeciesNames == MCMC.dat$Species[i]]
}
# tree is not utrametric (eventhough it is a BEAST output). Make that correction
# tree <- force.ultrametric(tree)
# remove unwanted tips
tree <- drop.tip(tree, processed.dat$SpecisName[is.na(processed.dat$male2N)])

# vector of chromosome number
chroms <- MCMC.dat$Chroms
names(chroms) <- MCMC.dat$Species

# arrange tip lables according to the taxa order of the tree
tip.labs <- matrix(data= NA,
                   nrow = Ntip(tree), 
                   ncol = 2)

for(i in 1:nrow(tip.labs)){
  tip.labs[i,1] <- tree$tip.label[i]
  tip.labs[i,2] <- MCMC.dat$HostType[MCMC.dat$Species == tree$tip.label[i]]
}

# make the contmap
map <- contMap(tree = tree,
        x = log(chroms),
        res = 100,
        lwd = 2,
        ftype = "off",
        type = "fan",outline = F, plot = F)

# change colours of contmap to viridis
map$cols <- viridis(n = 1001)
names(map$cols) <- 0:1000

# adgist the limits of the contMap so that it represents the real values
# not the log transformed values
map$lims <- c(min(chroms), max(chroms))

plot(map,
     type = "fan",
     res = 100,
     lwd = 3,
     outline = F,
     ftype = "off",
     fsize = c(0,.7),
     legend = F)

add.color.bar(50,
              map$cols,
              title="Haploid\nchromosome number",
              lims= map$lims,
              digits=3,
              prompt=FALSE,
              x=-100, # this is the possition where legend starts. We chose this
              # by using the locator function in R
              y=-91,
              lwd=4,
              fsize=.7,
              subtitle="")

# fix offset
.PlotPhyloEnv$last_plot.phylo$align.tip.label <- T

# plot tip lables
tiplabels(pch = 16, 
          col = c("blue","red")[as.numeric(tip.labs[,2]) + 1],
          cex = .5,
          offset = 1.5)

# plot legend for tip lables
points(x = rep(-100, 2),
       y = c(91,86),
       pch = 16, 
       cex = 1, 
       col = c("blue", "red"))

# add text
text(x = rep(-100, 2),
     y = c(91,86),
     labels = c("Generalist species", "Specialist species"),
     pos = 4,
     cex = .8)
