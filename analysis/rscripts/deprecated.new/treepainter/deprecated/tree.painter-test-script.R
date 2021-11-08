library(castor)
library(ape)
library(viridis)
library(phytools)
library(diversitree)

rate.class <- c(1,2,3,4,5)
rate <- c(.25,.5,1,2,4)
cols <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")


source("tree.painter.R")
source("get.rates.R")

tips <- 10

# make a dummy tree
# tree <- read.tree("../data/trees/papilionoidea-10.trees")
# tree <- tree[[1]]
# read in hosts
# hosts <- read.csv("../data/hosts/papilionoidea-hosts.csv")

#make tree
tree <- rcoal(n = tips)
tree$rates <- rep(3,Nedge(tree))
# make a simple model
qmat <- matrix(data = NA,
               nrow = 2,
               ncol = 2)
colnames(qmat) <- rownames(qmat) <- c("1","2")
qmat[1,1] <- -.4
qmat[1,2] <- .4
qmat[2,1] <- .4
qmat[2,2] <- -.4

# get tip states 
# trait <- sample(c(1,2), tips, replace = T)
trait <- sim.character(tree = tree,
                       pars = qmat,
                       model = "mkn",x0 = 2)

# make a simple model
qmat2 <- matrix(data = NA,
               nrow = 2,
               ncol = 2)
colnames(qmat2) <- rownames(qmat2) <- c("1","2")
qmat2[1,1] <- -.0018
qmat2[1,2] <- .0018
qmat2[2,1] <- .006
qmat2[2,2] <- -.006

# 
# painted.tree <- tree.paintR(tree = tree,
#                             tip_states = hosts$type + 1,
#                             qmat = qmat2,
#                             iter = 3,
#                             rate.class = rate.class,
#                             rate = rate,sample_mode = "median",return_br_table = T)

painted.tree <- tree.paintR.ver.3(tree = tree,
                  tip_states = trait,
                  qmat = qmat2,
                  iter = 10,
                  rate.class = rate.class,
                  rate = rate,
                  sample_mode = "median",return_br_table = T)

plot(tree,
     show.tip.label = F,
     edge.color = cols[as.factor(painted.tree$tree$rates)],
     edge.width = 2)

painted.tree$tree$rates

rate <- c(.25,.5,1,2,4)
get.rates(tree = painted.tree$tree, edge = 18, rate = rate)
