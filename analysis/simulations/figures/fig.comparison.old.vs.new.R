# Terrence Sylvester
# pradakshanas@gmail.com
# 22 February 2022

# read in data
dat.old.m <- read.csv("../results/power.old.method.csv", as.is = T)
dat.new.m <- read.csv("../results/power.new.method.csv", as.is = T)

# sort new method data by number of tips
dat.new.m <- dat.new.m[order(dat.new.m$ntip),]

# plot
par(mfrow = c(1,3))
for(i in c(1,2,5,10)){
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1),
       cex.axis = 1.5,
       cex.lab = 1.5)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 1)
  mtext(text = paste("Rates ratio =", i),side = 3,adj = 0,
        cex = 1)
  lines(x = dat.old.m$tips[dat.old.m$cond == 1 & dat.old.m$rr == i], y = dat.old.m$sigFission[dat.old.m$cond == 1 & dat.old.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 2 & dat.old.m$rr == i], y = dat.old.m$sigFission[dat.old.m$cond == 2 & dat.old.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 3 & dat.old.m$rr == i], y = dat.old.m$sigFission[dat.old.m$cond == 3 & dat.old.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 16, cex = 1.5)
  
  lines(x = dat.new.m$ntip[dat.new.m$cond == 1 & dat.new.m$rr == i], y = dat.new.m$fission[dat.new.m$cond == 1 & dat.new.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 2 & dat.new.m$rr == i], y = dat.new.m$fission[dat.new.m$cond == 2 & dat.new.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 3 & dat.new.m$rr == i], y = dat.new.m$fission[dat.new.m$cond == 3 & dat.new.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  
  # abline at power = 80%
  abline(h = 0.8, col  = "red", lty = 2)
  
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1),
       cex.axis = 1.5,
       cex.lab = 1.5)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 1)
  mtext(text = paste("Rates ratio =", i),side = 3,adj = 0,
        cex = 1)
  lines(x = dat.old.m$tips[dat.old.m$cond == 1 & dat.old.m$rr == i], y = dat.old.m$sigFusion[dat.old.m$cond == 1 & dat.old.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 2 & dat.old.m$rr == i], y = dat.old.m$sigFusion[dat.old.m$cond == 2 & dat.old.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 3 & dat.old.m$rr == i], y = dat.old.m$sigFusion[dat.old.m$cond == 3 & dat.old.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 16, cex = 1.5)
  
  lines(x = dat.new.m$ntip[dat.new.m$cond == 1 & dat.new.m$rr == i], y = dat.new.m$fusion[dat.new.m$cond == 1 & dat.new.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 2 & dat.new.m$rr == i], y = dat.new.m$fusion[dat.new.m$cond == 2 & dat.new.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 3 & dat.new.m$rr == i], y = dat.new.m$fusion[dat.new.m$cond == 3 & dat.new.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  
  # abline at power = 80%
  abline(h = 0.8, col  = "red", lty = 2)
  
  plot(x = NULL,
       y = NULL,
       xlab = "Number of tips",
       ylab = "Proportion of significant results",
       xlim = c(50,200),
       ylim = c(0,1),
       cex.axis = 1.5,
       cex.lab = 1.5)
  mtext(text = "Fission",side = 3,adj = 1,
        cex = 1)
  mtext(text = paste("Rates ratio =", i),side = 3,adj = 0,
        cex = 1)
  lines(x = dat.old.m$tips[dat.old.m$cond == 1 & dat.old.m$rr == i], y = dat.old.m$sigAneuploidy[dat.old.m$cond == 1 & dat.old.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 2 & dat.old.m$rr == i], y = dat.old.m$sigAneuploidy[dat.old.m$cond == 2 & dat.old.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 16, cex = 1.5)
  lines(x = dat.old.m$tips[dat.old.m$cond == 3 & dat.old.m$rr == i], y = dat.old.m$sigAneuploidy[dat.old.m$cond == 3 & dat.old.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 16, cex = 1.5)
  
  lines(x = dat.new.m$ntip[dat.new.m$cond == 1 & dat.new.m$rr == i], y = dat.new.m$aneuploidy[dat.new.m$cond == 1 & dat.new.m$rr == i],type = "b",col = "#e6610180", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 2 & dat.new.m$rr == i], y = dat.new.m$aneuploidy[dat.new.m$cond == 2 & dat.new.m$rr == i],type = "b",col = "#fdb86380", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  lines(x = dat.new.m$ntip[dat.new.m$cond == 3 & dat.new.m$rr == i], y = dat.new.m$aneuploidy[dat.new.m$cond == 3 & dat.new.m$rr == i],type = "b",col = "#5e3c9980", lwd = 2, pch = 17, lty = 2, cex = 1.5)
  
  # abline at power = 80%
  abline(h = 0.8, col  = "red", lty = 2)
  
  if(i == 10){
    # add legend
    legend("bottomright",inset=.02,
           legend=c("condition 1",
                    "condition 2",
                    "condition 3"),
           col=c("#e66101" ,"#fdb863" ,"#5e3c99"),
           lty=1,
           cex=1,
           lwd = 2,title = "Rates ratio")
  }
}


