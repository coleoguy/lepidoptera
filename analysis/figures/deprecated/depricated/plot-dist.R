asc.diff <- post.burnin$asc1 - post.burnin$asc2
desc.diff <- post.burnin$desc1 - post.burnin$desc2
pol.diff <- post.burnin$pol1 - post.burnin$pol2
tra.diff <- post.burnin$tran12 - post.burnin$tran21

plot(x = NULL, y = NULL, xlim = c(0,.4), ylim = c(-3,60), 
     xlab = expression(paste(Delta, "R"[x])), ylab = "Density")
polygon(density(asc.diff),col = rgb(0,0,1,.5), border = "Black", lwd = 2)
polygon(density(desc.diff),col = rgb(1,0,0,.5), border = "Black", lwd = 2)
polygon(density(pol.diff),col = rgb(0,1,0,.5), border = "Black", lwd = 2)

abline(v = 0, lty = 2, lwd = 2, col = "red")

points(x = rep(0.3,3), y = c(50,47,44), pch = 22, bg = c(rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(0,1,0,.5)), 
       cex = 2)
text(x = rep(0.3,3), y = c(50,47,44), 
     labels = c("Fusion", "Fission", "Polyploidy"), pos = 4)
a.hpd <- HPDinterval(as.mcmc(asc.diff))
d.hpd <- HPDinterval(as.mcmc(desc.diff))
p.hpd <- HPDinterval(as.mcmc(pol.diff))
t.hpd <- HPDinterval(as.mcmc(tra.diff))

segments(x0 = c(a.hpd[1], d.hpd[1], p.hpd[1]), 
         y0 = c(-1,-2,-3), 
         x1 = c(a.hpd[2], d.hpd[2], p.hpd[2]),
         y1 = c(-1,-2,-3),
         col = c(rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(0,1,0,.5)),
         lwd = 6)



plot(x = NULL, y = NULL, xlim = c(-.05,.02), ylim = c(-3,60), 
     xlab = expression(paste(Delta, "R"[x])), ylab = "Density")
polygon(density(tra.diff),col = rgb(0,0,1,.5), border = "Black", lwd = 2)
abline(v = 0, lty = 2, lwd = 2, col = "red")

segments(x0 = c(t.hpd[1]), 
         y0 = c(-3), 
         x1 = c(t.hpd[2]),
         y1 = c(-3),
         col = c(rgb(0,0,1,.5)),
         lwd = 6)






library(coda)
 
a.hpd <- HPDinterval(as.mcmc(asc.diff))
d.hpd <- HPDinterval(as.mcmc(desc.diff))
p.hpd <- HPDinterval(as.mcmc(pol.diff))

segments(x0 = c(a.hpd[1], d.hpd[1], p.hpd[1]), 
         y0 = c(-1,-2,-3), 
         x1 = c(a.hpd[2], d.hpd[2], p.hpd[2]),
         y1 = c(-1,-2,-3),
         col = c(rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(0,1,0,.5)),
         lwd = 6)


plot(x = NULL, y = NULL, xlim = c(0,.25), ylim = c(-.4,15), 
     xlab = expression(paste(Delta, "R")), ylab = "Density")
polygon(density(rnorm(n = 100, mean = .2,sd = .03)),col = rgb(1,0,0,.5), border = "red")
polygon(density(rnorm(n = 100, mean = .05,sd = .03)), col = rgb(0,0,1,.5), border = "blue")

points(x = rep(0.17,2), y = c(15,14), pch = 16, col = c("blue", "red"), 
       cex = 1.5)
text(x = rep(0.17,2), y = c(15,14), 
     labels = c("Generalist species", "Specialist species"), pos = 4)







expression(paste(delta, "R")







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
