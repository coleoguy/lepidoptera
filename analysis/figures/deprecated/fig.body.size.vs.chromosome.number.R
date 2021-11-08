# Terrence Sylvester
# pradakshanas@gmail.com
# 24 November 2020

# get helper functions
source("../rscripts/helper.functions.R")

# read in data
bsz <- read.csv("../data/body.size/body.size.csv", as.is = T)
bsz$speciesNames <- paste(bsz$Genus, bsz$Species, sep = "_")
# read in chromosome number data
dat <- read.delim("../data/chroms/chroms.txt", as.is = T)
# fill in species for data
dat$SpecisName <- paste(dat$other.names.genera,
                        dat$other.names.species,
                        sep = "_")
dat$SpecisName[dat$SpecisName == "_"] <- paste(dat$Genus[dat$SpecisName == "_"],
                                               dat$species[dat$SpecisName == "_"],
                                               sep = "_")
# sample a single chromosome number when there is a range
dat <- chromSampler(dat = dat)
# remove empty rows
bsz <- bsz[!is.na(bsz$wingspan.mm.),]
bsz <- bsz[!is.na(bsz$hostType),]
# remove duplicated species
bsz.uniqe.sp <- unique(bsz$speciesNames)
b.size <- matrix(data = NA,
                 nrow = length(bsz.uniqe.sp),
                 ncol = 3)
colnames(b.size) <- c("species", "wingspan", "hostType")
b.size <- as.data.frame(b.size)
for(i in 1:nrow(b.size)){
        b.size$species[i] <- bsz.uniqe.sp[i]
        b.size$wingspan[i] <- mean(bsz$wingspan.mm.[bsz$speciesNames == bsz.uniqe.sp[i]])
        b.size$hostType[i] <- sample(bsz$hostType[bsz$speciesNames == bsz.uniqe.sp[i]], 1)
}
# remove duplicated species
chrom.uniqe.sp <- unique(dat$SpecisName)
chroms <- matrix(data = NA,
                 nrow = length(chrom.uniqe.sp),
                 ncol = 2)
colnames(chroms) <- c("species", "hap")
chroms <- as.data.frame(chroms)
for(i in 1:nrow(chroms)){
        chroms$species[i] <- chrom.uniqe.sp[i]
        if(length(dat$male2N[dat$SpecisName == chrom.uniqe.sp[i]]) > 1){
                chroms$hap[i] <- sample(dat$male2N[dat$SpecisName == chrom.uniqe.sp[i]], 1)        
        }else{
                chroms$hap[i] <- dat$male2N[dat$SpecisName == chrom.uniqe.sp[i]]
        }
}
size <- b.size[b.size$species %in% chroms$species,]
chrom <- chroms[chroms$species %in% b.size$species,]
chrom$wingspan <- NA
chrom$host <- NA
for(i in 1:nrow(chrom)){
        chrom$wingspan[i] <- size$wingspan[size$species == chrom$species[i]]
        chrom$host[i] <- size$hostType[size$species == chrom$species[i]]
}
chrom$hap <- chrom$hap / 2
# remove all but chrom
rm(list = ls()[-4])
# plot
plot(y = chrom$hap,
     x = chrom$wingspan,
     col = c(rgb(1,0,0,.5), rgb(0,0,1,.5))[chrom$host + 1],
     pch = 16,
     ylab = "Haploid chromosome number",
     xlab = "Wingspan (mm)")
# labels
text(x = rep(110, 2),
     y = c(80,75),
     labels = c("Generalist", "Specialist"),
     col = c(rgb(1,0,0,1), rgb(0,0,1,1)),
     pos = 4)
# linear model
lmod <- lm(chrom$hap ~ chrom$wingspan)
lmod.sum <- summary(lmod)
abline(lmod, lwd = 2)
text(x = 100,
     y = 36,
     labels =  paste("p-value", round(lmod.sum$coefficients[2,4],3)),
     pos = 4)

