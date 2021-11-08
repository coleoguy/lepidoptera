load("../results/w.poly.allmatches.RData")
# evaluate burnin
plot(x[[1]]$p,type="l",ylim=c(-880,-760))
for(i in 2:100){
  lines(x[[i]]$p, col=rainbow(100)[i])
}
# looks like we can sample last 50%
results <- x[[1]][26:50,]
for(i in 2:100){
  results <- rbind(results, x[[i]][26:50,])
}
# evaluate impact of priors exp(rate=1)
# and upper limit upper=50
plot(density(results$asc1))
plot(density(results$asc2))
plot(density(results$desc1))
plot(density(results$desc2))
plot(density(results$asc1))
plot(density(results$pol1))
plot(density(results$pol2))

# convert to MY
