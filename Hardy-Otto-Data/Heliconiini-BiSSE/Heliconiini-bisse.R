library(diversitree)
library(geiger)

#load('musse-tree.Rda')
trees <- read.nexus('heliconiini-100.nex')
d <- read.csv('heliconiini-bisse.csv')
range <- d$range
names(range) <- d$sp

models.full <- list()
models.l <- list()

for (i in 1:100){
    print(i)
    phy <- trees[[i]]
    nc <- name.check(phy, range)
    t <- drop.tip(phy, nc$Tree.not.data)
    lik <- make.bisse(t, range)
    p <- starting.point.bisse(t)
    fit <- find.mle(lik, p)
    models.full <- c(models.full, list(fit))
    
    lik.l <- constrain(lik, lambda1 ~ lambda0)
    fit.l <- find.mle(lik.l, p[argnames(lik.l)])
    models.l <- c(models.l, list(fit.l))
    }

p.values <- vector()
for (i in 1:100){
    this.lrt <- anova(models.full[[i]], models.l[[i]])
    this.p <- this.lrt$"Pr(>|Chi|)"[[2]]
    p.values <- c(p.values, this.p)
    }
 
sumPump <- function(x, y){
    return(x$par[[y]])
    }
l0s <- lapply(models.full, sumPump, y=1)
mean(as.numeric(l0s))
sd(as.numeric(l0s))