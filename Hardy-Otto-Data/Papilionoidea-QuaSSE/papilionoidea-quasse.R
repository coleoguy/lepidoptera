library(diversitree)
library(geiger)

t <- read.tree('papilionoidea-MCC.tre')
d <- read.csv('papilionoidea-MCC.txt')
hosts <- d$hosts
names(hosts) <- d$sp

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)

set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=500, x0=1, single.lineage=FALSE)
states <- phy$tip.state
lik <- make.quasse(phy, states, sd(states), sigmoid.x, constant.x)
p <- starting.point.quasse(phy, states)

xr <- range(hosts) + c(-1,1) * 20 * p['diffusion']
linear.x <- make.linear.x(xr[1], xr[2])

make.butter <- function(lambda, mu){
    make.quasse(t, hosts, sd(hosts), lambda, mu, sampling.f=0.066)
    }
nodrift <- function(f){
    constrain(f, drift ~ 0)
    }

f.c <- make.butter(constant.x, constant.x)
f.l <- make.butter(linear.x, constant.x)

control <- list(parscale=.1, relto=0.001)
mle.c <- find.mle(nodrift(f.c), p, lower=0, control=control, verbose=0)
p.c <- mle.c$par
p.l <- c(p.c[1], l.m=0, p.c[2:3])

mle.l <- find.mle(nodrift(f.l), p.l, control=control, verbose=0)
