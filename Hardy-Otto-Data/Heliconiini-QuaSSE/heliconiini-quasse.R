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
    pardf <- data.frame(matrix(unlist(par.stack), nrow=100, byrow=T))
    colnames(pardf) <- names(par.stack[[1]])
    par.means <- colMeans(pardf)
    return(par.means)
    }

#get the trees and data
trees <- read.nexus('heliconiini-100.nex')
d <- read.csv('heliconiini-pds.txt')
pds <- d$pd
names(pds) <- d$sp

#a bit of ml search tuning
control <- list(parscale=.1, relto=0.001)

models.c <- list()
models.l <- list()
#loop
for (i in 1:100){
    print(i)
    raw.tree <- trees[[i]]
    td <- treedata(raw.tree, pds)
    t <- td$phy
    lik <- make.quasse(phy, pds, sd(pds), sigmoid.x, constant.x, control=control,sampling.f=0.29)
    p <- starting.point.quasse(t, pds)
    xr <- range(pds) + c(-1,1) * 20 * p['diffusion']
    linear.x <- make.linear.x(xr[1], xr[2])
    make.butter <- function(lambda, mu){
    make.quasse(t, pds, sd(pds), lambda, mu)
    }
    nodrift <- function(f){
    diversitree::constrain(f, drift ~ 0)
    }
    lik.c <- make.butter(constant.x, constant.x)
    lik.l <- make.butter(linear.x, constant.x)
    fit.c <- find.mle(nodrift(lik.c), p, lower=0, control=control, verbose=0)
    p.c <- fit.c$par
    p.l <- c(p.c[1], l.m=0, p.c[2:3])
    fit.l <- find.mle(nodrift(lik.l), p.l, control=control, verbose=0)
    models.c <- c(models.c, list(fit.c))
    models.l <- c(models.l, list(fit.l))
    }
    


p.values <- vector()
for (i in 1:100){
    this.lrt <- anova(models.l[[i]], models.c[[i]])
    this.p <- this.lrt$"Pr(>|Chi|)"[[2]]
    p.values <- c(p.values, this.p)
    }
    
slopes <- vector()
for (i in 1:100){
    s <- models.l[[i]]$par[2]
    slopes <- c(slopes, s)
    }
    
pars.l <- lapply(models.l, parStack)
pars.l.mean <- parMeans(pars.l)
ln.l <- lapply(models.l, likStack)

pars.c <- lapply(models.c, parStack)
pars.c.mean <- parMeans(pars.c)
ln.c <- lapply(models.c, likStack)

median.p <- median(as.numeric(p.values))
sig.ps <- p.values[p.values <= 0.05]
sig.p.freq <- length(sig.ps) / length(p.values)
    