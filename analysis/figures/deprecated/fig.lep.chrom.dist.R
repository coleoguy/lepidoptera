source("../rscripts/helper.functions - Copy.R")

library(ape)
library(ggplot2)
library(ggraptR)
library(phytools)
# read in data
dat <- read.delim("../data/chroms/chroms-working_copy.txt", as.is = T)
trees <- read.tree("../data/trees/processed.trees.new")
host <- read.csv("../data/hosts/papilionoidea-hosts.csv", as.is = T)

# fill in species for data
dat$SpecisName <- paste(dat$other.names.genera,
                        dat$other.names.species,
                        sep = "_")
dat$SpecisName[dat$SpecisName == "_"] <- paste(dat$Genus[dat$SpecisName == "_"],
                                               dat$species[dat$SpecisName == "_"],
                                               sep = "_")

# combine genus and species names to get the full name in hosts dataset and in
# chroms dataset
host$SpeciesNames <- paste(host$genus, host$species, sep = "_")
dat <- chromSampler(dat)
dat$male2N <- dat$male2N / 2

# remove NA
dat <- dat[!(is.na(dat$male2N)),]

ggplot(data = dat, aes(y = Family, x = male2N)) +
  geom_tile(stat="bin2d", position="identity", alpha=1, bins=200) +
  scale_x_continuous(breaks = seq(from = 0, to = 230, by = 10)) +
  theme_grey() + 
  scale_fill_viridis_c(trans = "log10", name = "Count")+
  xlab("Haploid chromosome number") + 
  ylab("Family")+
  theme(text=element_text(family="sans", face="plain", color="#000000", size=9, hjust=0.5, vjust=0.5))



##### some tests
# sample a single chromosome number when there is a range
new.dat <- chromSampler(dat = dat)
# remove species that lack male chromosome number data
new.dat <-  new.dat[!(is.na(new.dat$male2N)),]
# get the overlap between the traits and the phylogeny
overlap <- sp.matches(dat = new.dat, trees = trees)
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

MCMC.dat$Chroms <- MCMC.dat$Chroms * 2

# plot the distribution
# make the plot region
plot(x = NULL,
     y = NULL,
     xlim = c(0,1),
     ylim = c(0, max(MCMC.dat$Chroms)),
     xlab = "Host type",
     ylab = "Haploid chromosome count",xaxs = "i",
     xaxt = "n")

# plot tick marks
axis(side = 1,
     at= c(0.2,0.8), 
     labels= c("Generalist species", "Specialist species"))

# plot generalists
points(x = jitter(rep(.2, sum(MCMC.dat$HostType == 0)), 12),
       y = MCMC.dat$Chroms[MCMC.dat$HostType == 0],
       pch = 16,
       col = rgb(0,0,0,.5))

# plot specialists
points(x = jitter(rep(.8, sum(MCMC.dat$HostType == 1)), 3),
       y = MCMC.dat$Chroms[MCMC.dat$HostType == 1],
       pch = 16,
       col = rgb(0,0,0,.5))


var.test(MCMC.dat$Chroms ~ MCMC.dat$HostType)

var(MCMC.dat$Chroms[MCMC.dat$HostType == 1])
var(MCMC.dat$Chroms[MCMC.dat$HostType == 0])

chroms.1 <- (MCMC.dat$Chroms[MCMC.dat$HostType == 1])
names(chroms.1) <- MCMC.dat$Species[MCMC.dat$HostType == 1]

chroms.all <- MCMC.dat$Chroms
chroms.groups <- MCMC.dat$HostType

names(chroms.all) <- names(chroms.groups) <- MCMC.dat$Species

chroms.0 <- (MCMC.dat$Chroms[MCMC.dat$HostType == 0])
names(chroms.0) <- MCMC.dat$Species[MCMC.dat$HostType == 0]

t.test(MCMC.dat$Chroms ~ MCMC.dat$HostType)

phylANOVA(tree = tree,
          x = chroms.groups,
          y = chroms.all)

sum(is.na(MCMC.dat$Chroms[MCMC.dat$HostType == 0]))

tree <- force.ultrametric(tree)

contMap(tree = tree, x = chroms.all,, ftype = "off", lwd = 2,outline = F, type = "fan")
tiplabels(chroms.groups, offset = 3)
