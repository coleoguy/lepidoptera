# read in data
dat <- read.csv("../data/body.size/bsize-phylo-match.csv", as.is = T)
# remove data with no information for body size
dat <- dat[!is.na(dat$wingspanMean.mm.),]
# get the range
range(dat$wingspanMean.mm.)
# plot the distribution of mean body size
hist(dat$wingspanMean.mm.,breaks = 50)
# select the cutoff for high body size class and get the number of taxa for
# high body size class
sum(dat$wingspanMean.mm. > 60)
