# plotting

# load results
load("../results/all.matches.MCMC.RData")

par(mfcol = c(5,2), mar = c(3,3,3,3))

plot(results[[1]]$desc2, type = "l")

for(i in 2:10){
  plot(results[[i]]$desc2, type = "l", col = viridis::viridis(n = 10)[i])
}

post.burnin <- c()

for (i in 1:10) {
  x <- results[[i]][51:100,]
  post.burnin <- rbind(post.burnin, x)
}

plot(x = NULL, y = NULL, xlim = c(0,.25), ylim = c(0,100), 
     xlab = "Rate of fission (MYA)", ylab = "Density")
polygon(density(post.burnin$asc1),col = rgb(1,0,0,.5), border = "red")
polygon(density(post.burnin$asc2), col = rgb(0,0,1,.5), border = "blue")

points(x = rep(0.1875,2), y = c(100,90), pch = 16, col = c("blue", "red"), 
       cex = 1.5)
text(x = rep(0.1875,2), y = c(100,90), 
     labels = c("Generalist species", "Specialist species"), pos = 4)


plot(x = NULL, y = NULL, xlim = c(0,.4), ylim = c(0,200),
     xlab = "Rate of fusion (MYA)", ylab = "Density")
polygon(density(post.burnin$desc1),col = rgb(1,0,0,.5), border = "red")
polygon(density(post.burnin$desc2), col = rgb(0,0,1,.5), border = "blue")

points(x = rep(.3,2), y = c(200,180), pch = 16, col = c("blue", "red"), 
       cex = 1.5)
text(x = rep(.3,2), y = c(200,180), 
     labels = c("Generalist species", "Specialist species"), pos = 4)

plot(x = NULL, y = NULL, xlim = c(0,.4), ylim = c(0,200))
polygon(density(post.burnin$tran12),col = rgb(1,0,0,.5), border = "red")
polygon(density(post.burnin$tran21), col = rgb(0,0,1,.5), border = "blue")
