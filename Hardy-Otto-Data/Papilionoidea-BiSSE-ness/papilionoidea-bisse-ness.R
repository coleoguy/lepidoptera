library(diversitree)
library(geiger)

t <- read.tree('papilionoidea-MCC.tre')
d <- read.csv('papilionoidea-bisse-data.txt')
hosts <- d$hosts
names(hosts) <- d$sp

set.seed(3)

lik <- make.bisseness(t, hosts, sampling.f=c(0.083, 0.058))
startp <- starting.point.bisse(t)
bt <- branching.times(t)
tryq <- 1/2 * startp[['q01']] * sum(bt)/length(bt)
p <- c(startp[1:4], startp[5:6]/2, p0c=tryq, p0a=0.5, p1c=tryq, p1a=0.5)

lik.full <- constrain(lik, p0a ~ 0.5, p1a ~ 0.5)
lik.pc <- constrain(lik.full, p1c ~ p0c)
lik.pc0 <- constrain(lik, p0a ~ 0, p1a ~ 0, p0c ~ 0, p1c ~ 0)

fit <- find.mle(lik.full, p)
fit.pc <- find.mle(lik.pc, p[argnames(lik.pc)])
fit.pc0 <- find.mle(lik.pc0, p[argnames(lik.pc0)])