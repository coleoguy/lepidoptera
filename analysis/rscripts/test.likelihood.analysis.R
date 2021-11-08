# load library
library(chromePlus)
library(diversitree)

# load data
data("tree")
pars <- c(0.2, 0.2, 0, 0.2, 20)
# simulate chromosome numbers
chroms <-  simChrom(tree, pars=pars,
                    limits = c(1, 100), model = "2010") 
# make a data table with species names and chromosome number
dat <- as.data.frame(matrix(data = NA,nrow = length(chroms),ncol = 2))
colnames(dat) <- c("species", "chroms")
# fill the data table
dat$species <- names(chroms)
dat$chroms <- chroms

# get the chromosome matrix
chrom.mat <-  datatoMatrix(x = dat,
                           range = range(dat$chroms),
                           hyper = F)
# get the likelihood function
# make the likelihood function
lik <- make.mkn(tree = tree,
                states = chrom.mat,
                k = ncol(chrom.mat),
                strict = FALSE,
                control = list(method = "ode"))
# constrain the likelihood function
# constrain the likelihood function
con.lik <- constrainMkn(data = chrom.mat,
                        lik = lik,
                        polyploidy = F,
                        hyper = F,
                        # s.lambda = T,
                        # s.mu = T,
                        oneway = F,
                        constrain = list(drop.demi = T,
                                         drop.poly = F))
full.tree.lik <- con.lik(pars[c(1,2,4)])

# get likelihood after dropping a tip
tip.names <- tree$tip.label
# make a table to store likelihoods
lik.table <- as.data.frame(matrix(data = NA,
                                  nrow = length(tip.names),
                                  ncol = 3))
colnames(lik.table) <- c("tip.number", "tip.name", "likelihood")
# fill in rest of the datatable
lik.table$tip.number <- 1:length(tip.names)
lik.table$tip.name <- tip.names

# calculate the likelihood by dropping a tip
for(i in 1:nrow(lik.table)){
  # drop a single tip
  sub.tree <- drop.tip(tree, tip.names[i])
  # remove that data from the data table
  sub.MCMC.dat <- dat[dat$species != tip.names[i],]
  # make the likelihood function
  # get the range of chromosome number
  sub.rng <- c(range(dat$chroms, na.rm = T)[1] - 1,
               range(dat$chroms, na.rm = T)[2] + 1)
  # our states are coded such that specialists will be state 1 and
  # generalists will be state 2
  sub.chrom.mat <- datatoMatrix(x = sub.MCMC.dat,
                                range = sub.rng,
                                hyper = F)
  # make the likelihood function
  sub.lik <- make.mkn(tree = sub.tree,
                      states = sub.chrom.mat,
                      k = ncol(sub.chrom.mat),
                      strict = FALSE,
                      control = list(method = "ode"))
  # constrain the likelihood function
  sub.con.lik <- constrainMkn(data = sub.chrom.mat,
                              lik = sub.lik,
                              polyploidy = F,
                              hyper = F,
                              # s.lambda = T,
                              # s.mu = T,
                              oneway = F,
                              constrain = list(drop.demi = T,
                                               drop.poly = F))
  # get likelihood
  lik.table$likelihood[i] <-  sub.con.lik(pars[c(1,2,4)])
}
# plot the results
plot(lik.table$likelihood,
     cex = .5,
     pch = 16,
     xlab = "Dropped tip",
     ylab = "Likelihood",
     col = c("gray", "red")[(lik.table$likelihood > quantile(x = lik.table$likelihood, 0.95)) + 1])
# cutoff mark
abline(h = quantile(x = lik.table$likelihood, 0.95),
       col = "black",
       lty = 2)
text(x = 100,
     y = quantile(x = lik.table$likelihood, 0.95)-0.4,
     label = expression(paste("95"^"th","Quantile","")),
     pos = 3,
     cex = .7)

# get the likelihood of the unpruned tree
abline(h = full.tree.lik,
       col = "gray",
       lty = 2)
text(x = 95,
     y = full.tree.lik,
     label = "Likelihood of the unpruned tree",
     pos = 3,
     cex = .7)
#get the list of tips that have higher likelihood values
problamaticSPecies <- lik.table$tip.name[lik.table$likelihood > quantile(lik.table$likelihood, 0.95)]
# get the chromosome numers of these species
dat$chroms[dat$species %in% problamaticSPecies]
# see how the range changes
range(dat$chroms[!(dat$species %in% problamaticSPecies)])
range(dat$chroms)
