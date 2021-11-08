library(diversitree)
library(geiger)


parStack <- function(fit){
    this.par <- fit$par 
    return(this.par)
    }
    
likStack <- function(fit){
    this.lik <- fit$lnLik
    return(this.lik)
    }
    
parMeans <- function(par.stack){
    pardf <- data.frame(matrix(unlist(par.stack), nrow=10, byrow=T))
    colnames(pardf) <- names(par.stack[[1]])
    par.means <- colMeans(pardf)
    return(par.means)
    }



trees <- read.tree('papilionoidea-10.trees')
d <- read.csv('papilionoidea-bisse-data.txt')
range <- d$hosts
names(range) <- d$sp

models.full <- list()
models.l <- list()
models.m <- list()
models.r <- list()

for (i in 1:10){
    print(i)
    t <- trees[[i]]
    #phy <- trees[[i]]
    #td <- treedata(phy, range)
    #t <- td$phy
    lik <- make.bisse(t, range, sampling.f=c(0.083, 0.058))
    p <- starting.point.bisse(t)
    fit <- find.mle(lik, p)
    models.full <- c(models.full, list(fit))
    
    lik.l <- diversitree::constrain(lik, lambda1 ~ lambda0)
    fit.l <- find.mle(lik.l, p[argnames(lik.l)])
    models.l <- c(models.l, list(fit.l))
    
    lik.m <- diversitree::constrain(lik, mu1 ~ mu0)
    fit.m <- find.mle(lik.m, p[argnames(lik.m)])
    models.m <- c(models.m, list(fit.m))
    
    lik.r <- diversitree::constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0)
    fit.r <- find.mle(lik.r, p[argnames(lik.r)])
    models.r <- c(models.r, list(fit.r))
    }

p.values.l <- vector()
for (i in 1:10){
    this.lrt <- anova(models.full[[i]], models.l[[i]])
    this.p <- this.lrt$"Pr(>|Chi|)"[[2]]
    p.values.l <- c(p.values.l, this.p)
    }
    
p.values.m <- vector()
for (i in 1:10){
    this.lrt <- anova(models.full[[i]], models.m[[i]])
    this.p <- this.lrt$"Pr(>|Chi|)"[[2]]
    p.values.m <- c(p.values.m, this.p)
    }
    
p.values.r <- vector()
for (i in 1:10){
    this.lrt <- anova(models.full[[i]], models.r[[i]])
    this.p <- this.lrt$"Pr(>|Chi|)"[[2]]
    p.values.r <- c(p.values.r, this.p)
    }
    
pars.full <- lapply(models.full, parStack)
pars.mean <- parMeans(pars.full)
ln.full <- lapply(models.full, likStack)

pars.l <- lapply(models.l, parStack)
pars.l.mean <- parMeans(pars.l)
ln.l <- lapply(models.l, likStack)

pars.m <- lapply(models.m, parStack)
pars.m.mean <- parMeans(pars.m)
ln.m <- lapply(models.m, likStack)

pars.r <- lapply(models.r, parStack)
pars.r.mean <- parMeans(pars.r)
ln.r <- lapply(models.r, likStack)

median.p.l <- median(as.numeric(p.values.l))
sig.ps.l <- p.values.l[p.values.l <= 0.05]
sig.p.freq.l <- length(sig.ps.l) / length(p.values.l)

median.p.m <- median(as.numeric(p.values.m))
sig.ps.m <- p.values.m[p.values.m <= 0.05]
sig.p.freq.m <- length(sig.ps.m) / length(p.values.m)

median.p.r <- median(as.numeric(p.values.r))
sig.ps.r <- p.values.r[p.values.r <= 0.05]
sig.p.freq.r <- length(sig.ps.r) / length(p.values.r)
