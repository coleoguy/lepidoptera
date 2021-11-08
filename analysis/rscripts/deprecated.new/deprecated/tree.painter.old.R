library(diversitree)
library(phytools)
library(chromePlus)

# make a tree
tree <- rcoal(n = 20)
plot(tree, show.tip.label = F)
tree$edge.length <- tree$edge.length/ max(branching.times(tree))
# simulate traits
# binary trait
qmat.binary <- matrix(data = NA,
                      ncol = 2,
                      nrow = 2)

qmat.binary[1,1] <- -.5
qmat.binary[1,2] <- .5
qmat.binary[2,1] <- 1
qmat.binary[2,2] <- -1

colnames(qmat.binary) <- rownames(qmat.binary) <- c(0,1)
qmat.binary

tree$tip.states.binary <- sim.character(tree = tree, 
                                        pars = c(.3,.1), 
                                        model = "mk2",
                                        x0 = 0)

tiplabels(col = c("red", "blue")[tree$tip.states.binary+1], 
          pch =16,frame = "none")

tree$tip.states.discrete <- simChrom(tree = tree,
                                     limits = c(5,25),
                                     pars = c(0.01,1,0,0,15),
                                     model = "2010")
tiplabels(tree$tip.states.discrete, frame = "none", offset = 0.04)

dat <- as.data.frame(matrix(data = NA,
                            nrow = Ntip(tree),
                            ncol = 3))

colnames(dat) <- c("sp", "chrom", "bi")
dat$sp <- tree$tip.label
dat$chrom <- tree$tip.states.discrete
dat$bi <- tree$tip.states.binary

dat.to.m <- datatoMatrix(x = dat,
                         hyper = T)

musse.lik <- make.musse(tree = tree,
                        states = dat.to.m,
                        k = ncol(dat.to.m),
                        control = list("ODE"),
                        strict = F)

con.musse.lik <- constrainMuSSE(data = dat.to.m,
                                lik = musse.lik,
                                hyper = T,
                                polyploidy = F,
                                s.lambda = T,
                                s.mu = T,
                                constrain = list(drop.demi = T))

con.musse.lik(pars = c(0.01,1,0.01,1,0.3,0.1,1,1))


musse.pars <- find.mle(func = con.musse.lik,x.init = c(0.01,1,0.01,1,0.3,0.1,1,1))$par
con.musse.lik(musse.pars)

terminal.branches <- c()
for(i in 1:nrow(tree$edge)){
  if(tree$edge[i,2] %in% tree$edge[,1]){
    
  }else{
    terminal.branches <- c(terminal.branches, i)
  }
}

for(i in 1:length(terminal.branches)){
  for(j in 2:10){
    phy.sb <- tree
    phy.sb$edge.length[terminal.branches[i]] <- phy.sb$edge.length[terminal.branches[i]] * 20
    # phy.sb <- force.ultrametric(phy.sb, method = "extend")
    Musse.lik.sb <- make.musse(tree = phy.sb,
                                states = dat.to.m,
                                k = ncol(dat.to.m),
                                control = list("ODE"),
                                strict = F)
    
    if(bisse.lik.sb(pars) > bisse.lik(pars)){
      print(paste("found it")) 
    }
  }
}
