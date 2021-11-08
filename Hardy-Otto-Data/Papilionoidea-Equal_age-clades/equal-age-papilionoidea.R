library(ape)
library(phangorn)
library(geiger)
library(adephylo)


tree <- read.tree('mono-papilionoidea-MCC.tre')
tree$node.label <- NULL
btimes <- branching.times(tree)
youngsters <- btimes[btimes < 25]
targets <- vector()
index = 0
for (i in 1:length(youngsters)){
    n <- as.numeric(names(youngsters[i]))
    mom <- Ancestors(tree, n, type='parent')
    if (as.numeric(btimes[as.character(mom)]) > 25){
        targets[index] <- names(youngsters[i])
        index <- index + 1
        }
    }

winter <- function(node){
    tc <- extract.clade(tree, node)
    return(tc)
    }

buckets <- as.numeric(targets)
suckers <- lapply(buckets, winter)

loosers <- list()
for (i in 1:length(suckers)){
    if (length(suckers[[i]]$tip.label) > 10){
        loosers[i] <- suckers[i]
        }
    }

d <- read.csv('mono-papilionoidea-hosts.txt')
hosts <- d$hosts
names(hosts) <- d$sp

#get full tree
full.tree <- read.tree('papilionoidea-MCC.tre')
full.tree$node.label <- NULL

host_jumper <- function(xtree){
    tryCatch({fit <- fitDiscrete(xtree, hosts, model='ER')
    q <- fit$Trait1$q
    distances <- distTips(xtree)
    dist.matrix <- as.matrix(distances)
    dist.max <- which.max(dist.matrix)
    k <- arrayInd(dist.max, dim(dist.matrix))
    sp1 <- rownames(dist.matrix)[k[,1]]
    sp2 <- colnames(dist.matrix)[k[,2]]
    big.popa <- mrca(full.tree)[sp1, sp2]
    full.xtree <- extract.clade(full.tree, big.popa)
    rate <- rate.estimate(phy=full.xtree)
    dimer <- c(q, rate)
    return(dimer)
    }, error = function(e){
        NULL
    })
    }

randomization_test <- function(q, r, sample_size){
    cor.test <- cor(r, q, method='pearson')
    cor.values <- numeric(100)
    for (i in 1:100){
        index <- sample(1:sample_size, size=sample_size, replace=F)
        r.null <- r[index]
        cor.values[i] <- cor(r.null, q, method='pearson')
    }
    cor.values <- abs(cor.values)
    stat <- mean(cor.values >= cor.test)
    return(stat)
    }

winners <- loosers[sapply(loosers, function(x) length(x) > 0)]

answers <- lapply(winners, host_jumper)
daf <- data.frame(answers)
q <- as.numeric(daf[1,])
r <- as.numeric(daf[2,])
#plot(q,r)
#abline(lm(r ~ q))

test <- randomization_test(q, r, length(r))
print(test)
        


