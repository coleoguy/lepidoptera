# plot likelihood graphs for each condition and examine the output

# load data
load("mcmc.simulate.data.rr-5.RData")

# define fixed parameters
burnin <- .5
iter <- 100
burn <- -1:-(iter*burnin)

#### condition 1 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-110,-30),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.1[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-160,-80),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.1[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-260,-110),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.1[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-280,-140),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.1[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-320,-190),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.1[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.1[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.1[[1]][[i]][burn,])
}
post.burnin.condition.1.tree50.rr.5 <- y
write.csv(post.burnin.condition.1.tree50.rr.5, "post.burnin.condition.1.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.1[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.1[[2]][[i]][burn,])
}
post.burnin.condition.1.tree100.rr.5 <- y
write.csv(post.burnin.condition.1.tree100.rr.5, "post.burnin.condition.1.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.1[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.1[[3]][[i]][burn,])
}
post.burnin.condition.1.tree150.rr.5 <- y
write.csv(post.burnin.condition.1.tree150.rr.5, "post.burnin.condition.1.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.1[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.1[[4]][[i]][burn,])
}
post.burnin.condition.1.tree200.rr.5 <- y
write.csv(post.burnin.condition.1.tree200.rr.5, "post.burnin.condition.1.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.1[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.1[[5]][[i]][burn,])
}
post.burnin.condition.1.tree250.rr.5 <- y
write.csv(post.burnin.condition.1.tree250.rr.5, "post.burnin.condition.1.tree250.rr.5.csv",row.names = F)



#### condition 2 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-110,-30),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.2[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-170,-80),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.2[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-300,-110),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.2[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-300,-140),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.2[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-400,-190),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.2[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.2[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.2[[1]][[i]][burn,])
}
post.burnin.condition.2.tree50.rr.5 <- y
write.csv(post.burnin.condition.2.tree50.rr.5, "post.burnin.condition.2.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.2[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.2[[2]][[i]][burn,])
}
post.burnin.condition.2.tree100.rr.5 <- y
write.csv(post.burnin.condition.2.tree100.rr.5, "post.burnin.condition.2.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.2[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.2[[3]][[i]][burn,])
}
post.burnin.condition.2.tree150.rr.5 <- y
write.csv(post.burnin.condition.2.tree150.rr.5, "post.burnin.condition.2.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.2[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.2[[4]][[i]][burn,])
}
post.burnin.condition.2.tree200.rr.5 <- y
write.csv(post.burnin.condition.2.tree200.rr.5, "post.burnin.condition.2.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.2[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.2[[5]][[i]][burn,])
}
post.burnin.condition.2.tree250.rr.5 <- y
write.csv(post.burnin.condition.2.tree250.rr.5, "post.burnin.condition.2.tree250.rr.5.csv",row.names = F)

#### condition 3 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-110,-30),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.3[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-170,-50),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.3[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-250,-100),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.3[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-300,-140),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.3[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-400,-190),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.3[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.3[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.3[[1]][[i]][burn,])
}
post.burnin.condition.3.tree50.rr.5 <- y
write.csv(post.burnin.condition.3.tree50.rr.5, "post.burnin.condition.3.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.3[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.3[[2]][[i]][burn,])
}
post.burnin.condition.3.tree100.rr.5 <- y
write.csv(post.burnin.condition.3.tree100.rr.5, "post.burnin.condition.3.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.3[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.3[[3]][[i]][burn,])
}
post.burnin.condition.3.tree150.rr.5 <- y
write.csv(post.burnin.condition.3.tree150.rr.5, "post.burnin.condition.3.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.3[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.3[[4]][[i]][burn,])
}
post.burnin.condition.3.tree200.rr.5 <- y
write.csv(post.burnin.condition.3.tree200.rr.5, "post.burnin.condition.3.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.3[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.3[[5]][[i]][burn,])
}
post.burnin.condition.3.tree250.rr.5 <- y
write.csv(post.burnin.condition.3.tree250.rr.5, "post.burnin.condition.3.tree250.rr.5.csv",row.names = F)

#### condition 4 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-110,-30),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.4[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-200,-80),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.4[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-300,-110),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.4[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-300,-140),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.4[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-400,-190),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.4[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.4[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.4[[1]][[i]][burn,])
}
post.burnin.condition.4.tree50.rr.5 <- y
write.csv(post.burnin.condition.4.tree50.rr.5, "post.burnin.condition.4.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.4[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.4[[2]][[i]][burn,])
}
post.burnin.condition.4.tree100.rr.5 <- y
write.csv(post.burnin.condition.4.tree100.rr.5, "post.burnin.condition.4.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.4[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.4[[3]][[i]][burn,])
}
post.burnin.condition.4.tree150.rr.5 <- y
write.csv(post.burnin.condition.4.tree150.rr.5, "post.burnin.condition.4.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.4[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.4[[4]][[i]][burn,])
}
post.burnin.condition.4.tree200.rr.5 <- y
write.csv(post.burnin.condition.4.tree200.rr.5, "post.burnin.condition.4.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.4[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.4[[5]][[i]][burn,])
}
post.burnin.condition.4.tree250.rr.5 <- y
write.csv(post.burnin.condition.4.tree250.rr.5, "post.burnin.condition.4.tree250.rr.5.csv",row.names = F)

#### condition 5 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.5[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.5[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.5[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.5[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.5[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.5[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.5[[1]][[i]][burn,])
}
post.burnin.condition.5.tree50.rr.5 <- y
write.csv(post.burnin.condition.5.tree50.rr.5, "post.burnin.condition.5.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.5[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.5[[2]][[i]][burn,])
}
post.burnin.condition.5.tree100.rr.5 <- y
write.csv(post.burnin.condition.5.tree100.rr.5, "post.burnin.condition.5.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.5[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.5[[3]][[i]][burn,])
}
post.burnin.condition.5.tree150.rr.5 <- y
write.csv(post.burnin.condition.5.tree150.rr.5, "post.burnin.condition.5.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.5[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.5[[4]][[i]][burn,])
}
post.burnin.condition.5.tree200.rr.5 <- y
write.csv(post.burnin.condition.5.tree200.rr.5, "post.burnin.condition.5.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.5[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.5[[5]][[i]][burn,])
}
post.burnin.condition.5.tree250.rr.5 <- y
write.csv(post.burnin.condition.5.tree250.rr.5, "post.burnin.condition.5.tree250.rr.5.csv",row.names = F)

#### condition 6 ####
#### plot likelihood ####
# 50 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.6[[1]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 100 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.6[[2]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 150 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.6[[3]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 200 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.6[[4]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

# 250 tips
plot(x = NULL,
     y = NULL,
     xlim = c(0,100),
     ylim = c(-350,-150),
     xlab = "Generation",
     ylab = "Likelihood")

for(i in 1:100){
  lines(results.condition.6[[5]][[i]][[8]], col = rainbow(100, alpha = 0.5)[i])  
}

#### get burnin for each tree ####
#50 tips
y <- c()
y <- results.condition.6[[1]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.6[[1]][[i]][burn,])
}
post.burnin.condition.6.tree50.rr.5 <- y
write.csv(post.burnin.condition.6.tree50.rr.5, "post.burnin.condition.6.tree50.rr.5.csv",row.names = F)

#100 tips
y <- c()
y <- results.condition.6[[2]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.6[[2]][[i]][burn,])
}
post.burnin.condition.6.tree100.rr.5 <- y
write.csv(post.burnin.condition.6.tree100.rr.5, "post.burnin.condition.6.tree100.rr.5.csv",row.names = F)

#150 tips
y <- c()
y <- results.condition.6[[3]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.6[[3]][[i]][burn,])
}
post.burnin.condition.6.tree150.rr.5 <- y
write.csv(post.burnin.condition.6.tree150.rr.5, "post.burnin.condition.6.tree150.rr.5.csv",row.names = F)

#200 tips
y <- c()
y <- results.condition.6[[4]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.6[[4]][[i]][burn,])
}
post.burnin.condition.6.tree200.rr.5 <- y
write.csv(post.burnin.condition.6.tree200.rr.5, "post.burnin.condition.6.tree200.rr.5.csv",row.names = F)

#250 tips
y <- c()
y <- results.condition.6[[5]][[1]][burn,]
for(i in 2:iter){
  y <- rbind(y, results.condition.6[[5]][[i]][burn,])
}
post.burnin.condition.6.tree250.rr.5 <- y
write.csv(post.burnin.condition.6.tree250.rr.5, "post.burnin.condition.6.tree250.rr.5.csv",row.names = F)

