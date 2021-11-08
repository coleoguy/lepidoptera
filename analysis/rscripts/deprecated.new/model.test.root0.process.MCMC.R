# Terrence Sylvester
# 28th October 2020
# pradakshanas@gmail.com
# process and visualise the mcmc output

# load results
load("../results/model.test.root.0.RData")
results <- x

# remove unwanted objects
ls()

# lets fix iteration and liklikood values
for(i in 1:length(results)){
  multiplier <- 1/results[[i]]$i[1]
  results[[i]]$i <- results[[i]]$i * multiplier
  results[[i]]$p <- results[[i]]$p * multiplier
}

# see when the mcmc have converged
plot(results[[1]]$p, type = "l",
     ylim = c(-1500,-1050),
     col = rainbow(1,alpha = .5))

for(i in 2:100){
  lines(results[[i]]$p, col = rainbow(n = 100,alpha = .5)[i])
}

# looks like mcmc have converged around 20 generations

par(mfrow = c(2,3))
for(i in 1:100){
  plot(results[[i]]$asc1,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "specialists gains")
  
  plot(results[[i]]$desc1,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "specialists loss")
  
  plot(results[[i]]$pol1,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "specialists poly")
  
  plot(results[[i]]$asc2,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "generalists gains")
  
  plot(results[[i]]$desc2,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "generalists loss")
  
  plot(results[[i]]$pol2,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "generalists poly")
  
}

# all trees look good

# see when the mcmc have converged
plot(results[[1]]$p, type = "l",
     ylim = c(-1700,-1450),
     col = rainbow(1,alpha = .5))

par(mfrow = c(1,10))
for(i in 1:100){
  plot(results[[i]]$p,
       # col = rainbow(n = 100,alpha = .5)[i], 
       type = "l",
       main = paste("iteration", i),
       sub = paste("Tree", tree.num[i]),
       ylab = "",
       xlab = "")
}

# Inspection of the likelihood values also shows that in above mentioned trees
# the liklihood is still increasing.

# lets discard these trees from our resutls
proc.results <- results

# sample the burnin portion
# burnin 
burnin <- .50
iter <- 100
burn <- -1:-(iter*burnin)
y <- proc.results[[1]][burn,]
for (i in c(2,c(4:length(proc.results)))) {
  y <- rbind(y, proc.results[[i]][burn,])
}
post.burnin <- y

# remove all unwanted results
rm(list = ls()[-10])

# save results
write.csv(x = post.burnin,
          file = "../results/model.test.root.0.post.burnin.csv",
          row.names = F)
